import pickle, numpy as np, re
from pathlib import Path
from config import RUNS, WORK
P3 = WORK
D = pickle.load(open(P3 / 'b11_res.pkl', 'rb'))
RUN = RUNS


def mon(case):
    t = []; m = []
    for l in open(RUN / case / 'postProcessing/maxSpeed/0/fieldMinMax.dat'):
        if l.startswith('#') or not l.strip(): continue
        g = re.match(r'\s*([0-9.e+-]+)\s+mag\(U\)\s+([0-9.e+-]+)\s+\([^)]*\)\s+([0-9.e+-]+)\s+\(', l)
        t.append(float(g[1])); m.append(float(g[3]))
    return np.array(t), np.array(m)


a, sa = D['V1a_static']; b, sb = D['S_rhoST']


def col(d, q, sel=None):
    v = np.array([np.nan if r[q] is None else r[q] for r in d['snaps']])
    return v if sel is None else v[sel]


T = np.array([r['t'] for r in a['snaps']])
w1 = T <= 1.0 + 1e-3; w03 = T >= 0.3; wall = T > 0


def rg(v): return f"{np.nanmin(v):.3f} to {np.nanmax(v):.3f}"
def rt(x, y): return f"{x / y:.3f}"


L = []
add = L.append
ta, ma = mon('V1a_static'); tb, mb = mon('S_rhoST')
k1 = (ta >= 0.3) & (ta <= 1.0001)
add("# B11: V1a_static (standard interFoam) against S_rhoST (density-weighted surface tension)\n")
add("Numbers only. All values are computed by `p3/b11_srho.py` (log `p3/b11_run.log`, per-snapshot table `p3/b11_snap.csv`, report generator `p3/b11_report.py`). Cell values come from the written fields of `<P2_RUNS>/<run>/<t>/{alpha.water,U,p}` and `0/{C,V}`, read with `scal`/`vec` and `tb_post` exactly as `snap.py` does (liquid = alpha > 0.99, gas = alpha < 0.01, interface = the rest; `tb_post.effective_radius`, `tb_post.laplace_jump(Cx, Cy, alpha, p, 2.0, 0.25)` with the default margin 0.3; sigma = 0.0999 from transportProperties in both runs). maxSpeed monitor = column 'max' of `postProcessing/maxSpeed/0/fieldMinMax.dat`; liquidVolume = `postProcessing/liquidVolume/0/volFieldValue.dat`.\n")
add("Case check (diff of the two case folders, whitespace and line endings ignored): the only differences are `application interFoamRhoST` in system/controlDict and `rhoWeightedSurfaceTension yes;` in the PIMPLE dict of system/fvSolution; constant/ is identical, 0/C and 0/alpha.water are byte-identical. Both ran 23 061 steps to 3.95482 s with the same time step (monitor spacing 1.7142e-3 to 1.7151e-3 s per 10 steps), so the arithmetic mean of monitor samples equals the time mean (the trapezoid mean agrees to 3 digits, listed below). Mesh: 96 x 96 = 9216 cells, x uniform (dx = 0.04167), r stretched (dr = 0.0061 at r = 0.25), cell aspect ratio at the interface about 6.8. Boundaries: x = 0 wall, x = 4 open (totalPressure 0, pressureInletOutletVelocity), tube wall no-slip.\n")
# ---------- 1
add("## 1. Reproduction of the paper row for V1a_static (to t = 1 s)\n")
add("| Quantity | Paper (Table, stretched row) | My value, V1a_static | Source |\n|---|---|---|---|")
add(f"| Mean maximum speed, 0.3 to 1 s | 7.27 m/s | **{sa['ms_1s']:.3f} m/s** (n = {sa['ms_1s_n']} samples; trapezoid time mean {sa['ms_1s_trap']:.3f}) | maxSpeed monitor, samples from t = {ta[k1][0]:.4f} to {ta[k1][-1]:.4f} s |")
add(f"| Liquid peak, snapshots to 1 s | 0.41 to 0.60 m/s | **{sa['maxliq_1s']['min']:.3f} to {sa['maxliq_1s']['max']:.3f} m/s** (5 snapshots, t = 0.198 to 0.989 s) | max of abs(U) over cells with alpha > 0.99 |")
add(f"| Pressure jump, snapshots to 1 s | 0.72 to 0.78 Pa (Laplace 0.799) | **{sa['jump_1s']['min']:.3f} to {sa['jump_1s']['max']:.3f} Pa** (Laplace 2 sigma / Reff = {sa['target_1s']['mean']:.4f} Pa) | tb_post.laplace_jump |")
add("\nMatch: 7.268 rounds to 7.27; 0.408 to 0.605 rounds to 0.41 to 0.60; 0.7150 to 0.7822 rounds to 0.72 to 0.78 (0.71501 is just above the 0.715 rounding edge, so this end of the range is sensitive to the last digit). The method reproduces the paper row, so the same code was used on S_rhoST below.\n")
# ---------- 2
add("## 2. Same three quantities, same window (to t = 1 s, mean over 0.3 to 1 s)\n")
add("| Quantity | V1a_static | S_rhoST | S_rhoST / V1a_static |\n|---|---|---|---|")
add(f"| Mean maximum speed, 0.3 to 1 s (m/s) | {sa['ms_1s']:.3f} | {sb['ms_1s']:.3f} | {rt(sb['ms_1s'], sa['ms_1s'])} |")
add(f"| Max of the maximum speed in that window (m/s) | {sa['ms_1s_max']:.3f} | {sb['ms_1s_max']:.3f} | {rt(sb['ms_1s_max'], sa['ms_1s_max'])} |")
add(f"| Liquid peak, 5 snapshots (m/s) | {rg(col(a,'maxliq',w1))} (mean {sa['maxliq_1s']['mean']:.3f}) | {rg(col(b,'maxliq',w1))} (mean {sb['maxliq_1s']['mean']:.3f}) | mean {rt(sb['maxliq_1s']['mean'], sa['maxliq_1s']['mean'])}; min {rt(sb['maxliq_1s']['min'], sa['maxliq_1s']['min'])}; max {rt(sb['maxliq_1s']['max'], sa['maxliq_1s']['max'])} |")
add(f"| Gas peak, 5 snapshots (m/s) | {rg(col(a,'maxgas',w1))} (mean {sa['maxgas_1s']['mean']:.3f}) | {rg(col(b,'maxgas',w1))} (mean {sb['maxgas_1s']['mean']:.3f}) | mean {rt(sb['maxgas_1s']['mean'], sa['maxgas_1s']['mean'])} |")
add(f"| Pressure jump, 5 snapshots (Pa) | {rg(col(a,'jump',w1))} (mean {sa['jump_1s']['mean']:.3f}, {sa['jump_err_1s']:+.1f}% of Laplace) | {rg(col(b,'jump',w1))} (mean {sb['jump_1s']['mean']:.3f}, {sb['jump_err_1s']:+.1f}% of Laplace) | mean {rt(sb['jump_1s']['mean'], sa['jump_1s']['mean'])}; min {rt(sb['jump_1s']['min'], sa['jump_1s']['min'])}; max {rt(sb['jump_1s']['max'], sa['jump_1s']['max'])} |")
# ---------- 3
add("\n## 3. Whole run\n")
add(f"maxSpeed monitor: t = 0.3 to 3.9546 s ({sa['ms_all_n']} samples each). Snapshots: the 19 written snapshots with t >= 0.3 s (0.395 to 3.955 s) are the main rows, and the 20 snapshots from 0.198 s (the set behind the paper's 0.728 Pa) are given in a second block. The pressure jump of S_rhoST is undefined (the function returns None: no gas cell inside 0.7 R of (x = 2, r = 0)) in the 7 snapshots from t = 2.768 s on, because the bubble has moved away from x = 2 (see item 5 and the Expert note). Its jump statistics therefore use 12 snapshots (t >= 0.3) or 13 (all 20), against 19 or 20 for V1a. A like-for-like row on the same 12 snapshots is added.\n")
add("| Quantity | V1a_static | S_rhoST | S_rhoST / V1a_static |\n|---|---|---|---|")
add(f"| maxSpeed, time mean 0.3 to 3.95 s (m/s) | {sa['ms_all']:.3f} (trapezoid {sa['ms_all_trap']:.3f}) | {sb['ms_all']:.3f} (trapezoid {sb['ms_all_trap']:.3f}) | {rt(sb['ms_all'], sa['ms_all'])} |")
add(f"| maxSpeed, max over 0.3 to 3.95 s (m/s) | {sa['ms_all_max']:.3f} | {sb['ms_all_max']:.3f} | {rt(sb['ms_all_max'], sa['ms_all_max'])} |")
add(f"| maxSpeed, max over the whole run incl. start-up (m/s, at t) | {sa['ms_peak_run']:.3f} (t = {sa['t_peak_run']:.4f}) | {sb['ms_peak_run']:.3f} (t = {sb['t_peak_run']:.4f}) | {rt(sb['ms_peak_run'], sa['ms_peak_run'])} |")


def win(s): return '; '.join(f"{x0:g}-{x1:g} s: {v:.3f}" for x0, x1, v in s['mon_win'])


add(f"| maxSpeed, time mean per window (m/s) | {win(sa)} | {win(sb)} | {', '.join(rt(vb[2], va[2]) for va, vb in zip(sa['mon_win'], sb['mon_win']))} |")
for name, lab in [('ge03', 'snapshots t >= 0.3 s, n = 19'), ('all20', 'all 20 snapshots from 0.198 s')]:
    qa, qb = sa['maxliq_' + name], sb['maxliq_' + name]
    ga, gb = sa['maxgas_' + name], sb['maxgas_' + name]
    add(f"| Liquid peak, mean (min to max), {lab} (m/s) | {qa['mean']:.3f} ({qa['min']:.3f} to {qa['max']:.3f}) | {qb['mean']:.3f} ({qb['min']:.3f} to {qb['max']:.3f}) | mean {rt(qb['mean'], qa['mean'])}; min {rt(qb['min'], qa['min'])}; max {rt(qb['max'], qa['max'])} |")
    add(f"| Gas peak, mean (min to max), {lab} (m/s) | {ga['mean']:.3f} ({ga['min']:.3f} to {ga['max']:.3f}) | {gb['mean']:.3f} ({gb['min']:.3f} to {gb['max']:.3f}) | mean {rt(gb['mean'], ga['mean'])}; min {rt(gb['min'], ga['min'])}; max {rt(gb['max'], ga['max'])} |")
    ja, jb = sa['jump_' + name], sb['jump_' + name]
    add(f"| Pressure jump, mean (min to max), fixed centre x = 2, {lab} (Pa) | {ja['mean']:.3f} ({ja['min']:.3f} to {ja['max']:.3f}), n = {sa['jump_nvalid_' + name]}, {sa['jump_err_' + name]:+.1f}% of Laplace | {jb['mean']:.3f} ({jb['min']:.3f} to {jb['max']:.3f}), n = {sb['jump_nvalid_' + name]}, {sb['jump_err_' + name]:+.1f}% of Laplace | mean {rt(jb['mean'], ja['mean'])}; min {rt(jb['min'], ja['min'])}; max {rt(jb['max'], ja['max'])} |")
okb = ~np.isnan(col(b, 'jump')) & w03
ja_l = col(a, 'jump')[okb]; jb_l = col(b, 'jump')[okb]
add(f"| Pressure jump on the same {int(okb.sum())} snapshots (t = {T[okb][0]:.3f} to {T[okb][-1]:.3f} s) where S_rhoST is defined (Pa) | {ja_l.mean():.3f} ({ja_l.min():.3f} to {ja_l.max():.3f}), sd {ja_l.std(ddof=1):.3f} | {jb_l.mean():.3f} ({jb_l.min():.3f} to {jb_l.max():.3f}), sd {jb_l.std(ddof=1):.3f} | mean {rt(jb_l.mean(), ja_l.mean())}; min {rt(jb_l.min(), ja_l.min())}; max {rt(jb_l.max(), ja_l.max())} |")
ja_c, jb_c = sa['jumpc_ge03'], sb['jumpc_ge03']
add(f"| Supplement, NOT the paper definition: same function with the centre moved to the gas centroid, t >= 0.3 s, n = 19 (Pa) | {ja_c['mean']:.3f} ({ja_c['min']:.3f} to {ja_c['max']:.3f}), {ja_c['err']:+.1f}% | {jb_c['mean']:.3f} ({jb_c['min']:.3f} to {jb_c['max']:.3f}), {jb_c['err']:+.1f}% | mean {rt(jb_c['mean'], ja_c['mean'])} |")
add(f"| Laplace value 2 sigma / Reff, mean over snapshots t >= 0.3 s (Pa) | {sa['target_ge03']['mean']:.4f} | {sb['target_ge03']['mean']:.4f} | |")
# ---------- 4
add("\n## 4. Where the largest speed sits in the last snapshot (t = 3.95482 s)\n")
la, lb = sa['last'], sb['last']
add("| Quantity | V1a_static | S_rhoST |\n|---|---|---|")
add(f"| Largest speed (m/s) | {la['max']:.3f} | {lb['max']:.3f} |")
add(f"| Class of that cell | {la['cls']}, alpha = {la['alpha']:.2e} | {lb['cls']}, alpha = {lb['alpha']:.2e} |")
add(f"| Position x, r | {la['x']:.4f}, {la['r']:.4f} | {lb['x']:.4f}, {lb['r']:.4f} |")
add(f"| Velocity components (Ux, Ur) (m/s) | ({la['Ux']:.3f}, {la['Ur']:.3f}) | ({lb['Ux']:.3f}, {lb['Ur']:.3f}) |")
add(f"| Largest speed in gas / interface / liquid (m/s) | {sa['last_snap']['maxgas']:.3f} / {sa['last_snap']['maxitf']:.3f} / {sa['last_snap']['maxliq']:.3f} | {sb['last_snap']['maxgas']:.3f} / {sb['last_snap']['maxitf']:.3f} / {sb['last_snap']['maxliq']:.3f} |")
add(f"| Mean abs(U) over interface cells, unweighted (m/s) (n cells) | {sa['last_snap']['meanitf']:.3f} (n = {sa['last_snap']['nitf']}) | {sb['last_snap']['meanitf']:.3f} (n = {sb['last_snap']['nitf']}) |")
add(f"| Same, volume-weighted (m/s) | {sa['last_snap']['meanitf_vw']:.3f} | {sb['last_snap']['meanitf_vw']:.3f} |")
add(f"| Mean abs(U) over interface cells, mean of the 20 snapshots (m/s); range of cell count | {sa['meanitf_all']:.3f}; n = {sa['nitf_range'][0]} to {sa['nitf_range'][1]} | {sb['meanitf_all']:.3f}; n = {sb['nitf_range'][0]} to {sb['nitf_range'][1]} |")
add(f"| Monitor at its last sample (t = {sa['mon_t_last']:.4f} s): max speed at (x, r) | {sa['mon_last']:.3f} at ({sa['mon_last_loc'][0]:.4f}, {sa['mon_last_loc'][1]:.4f}) | {sb['mon_last']:.3f} at ({sb['mon_last_loc'][0]:.4f}, {sb['mon_last_loc'][1]:.4f}) |")
add("\nThe maximum is in the gas, next to the interface (r about 0.19 to 0.22 against a bubble radius of 0.25), in both runs. The monitor values differ slightly from the snapshot values only because the last monitor sample is 1.7e-4 s earlier than the last write.\n")
# ---------- 5
add("## 5. Shape\n")
ha, hb = sa['shape'], sb['shape']
sm_a = np.array([r['smeared'] for r in a['snaps']]); sm_b = np.array([r['smeared'] for r in b['snaps']])
add("| Quantity | V1a_static | S_rhoST |\n|---|---|---|")
add(f"| Reff at t = 0 (from 0/alpha.water) (m) | {ha['Reff0']:.6f} | {hb['Reff0']:.6f} |")
add(f"| Reff at first written snapshot, t = 0.1977 s (m) | {ha['Reff_first']:.6f} | {hb['Reff_first']:.6f} |")
add(f"| Reff at end, t = 3.9548 s (m) | {ha['Reff_last']:.6f} | {hb['Reff_last']:.6f} |")
add(f"| Reff, min to max over the 20 snapshots (m) | {ha['Reff_min']:.6f} to {ha['Reff_max']:.6f} | {hb['Reff_min']:.6f} to {hb['Reff_max']:.6f} |")
add(f"| Smeared gas share at t = 0 (%) | {float(ha['sm0']):.3f} | {float(hb['sm0']):.3f} |")
add(f"| Smeared gas share at t = 0.1977 s (%) | {float(ha['sm_first']):.3f} | {float(hb['sm_first']):.3f} |")
add(f"| Smeared gas share at end, t = 3.9548 s (%) | {float(ha['sm_last']):.3f} | {float(hb['sm_last']):.3f} |")
add(f"| Smeared gas share, mean (min to max) over 20 snapshots (%) | {sm_a.mean():.3f} ({sm_a.min():.3f} to {sm_a.max():.3f}) | {sm_b.mean():.3f} ({sm_b.min():.3f} to {sm_b.max():.3f}) |")
add(f"| Smeared share, mean of the last 5 snapshots, t = 3.17 to 3.95 s (%) | {sm_a[-5:].mean():.3f} | {sm_b[-5:].mean():.3f} |")
add(f"| Gas centroid x at end (m); start = 2.0000 | {sa['xc_last']:.4f} | {sb['xc_last']:.4f} |")
add(f"| Gas centroid x, min to max over the run (m) | {sa['xc_range'][0]:.4f} to {sa['xc_range'][1]:.4f} | {sb['xc_range'][0]:.4f} to {sb['xc_range'][1]:.4f} |")
add(f"| Centroid velocity, linear fit over snapshots t >= 2.3 s (m/s) | {sa['xc_slope_ge23']:+.3f} | {sb['xc_slope_ge23']:+.3f} |")
add(f"| Extent of cells with alpha < 0.5 in x, min to max over snapshots (m) | {sa['ext_x_range5'][0]:.3f} to {sa['ext_x_range5'][1]:.3f} | {sb['ext_x_range5'][0]:.3f} to {sb['ext_x_range5'][1]:.3f} |")
add(f"| Largest r of cells with alpha < 0.5, min to max over snapshots (m) | {sa['ext_r_range5'][0]:.3f} to {sa['ext_r_range5'][1]:.3f} | {sb['ext_r_range5'][0]:.3f} to {sb['ext_r_range5'][1]:.3f} |")
add("\nReff is computed from the gas volume (tb_post.effective_radius), so it only tracks volume (its snapshot-to-snapshot change is the gas-volume drift of item 6) and says nothing about shape. The smeared share starts at 0 because setFields marks whole cells. Centroid and extent rows come from `b11_srho.py` (the extent of alpha < 0.5 uses cell centres, so values are quantised to the cell size, 0.0417 in x and about 0.006 in r); `b11_drift.py` and `b11_ux.py` (logs in this folder if kept) were used for the cross-checks quoted in the Expert note.\n")
# ---------- 6
add("## 6. Gas volume drift (wedge volume, m^3)\n")
add("| Quantity | V1a_static | S_rhoST |\n|---|---|---|")
add(f"| Gas volume at t = 0, snapshot 0/alpha.water | {sa['gas0_snap']:.8e} | {sb['gas0_snap']:.8e} |")
add(f"| Gas volume at 3.95482 s, snapshot | {sa['gasV_last']:.8e} | {sb['gasV_last']:.8e} |")
add(f"| Drift, snapshot last vs t = 0 (%) | {sa['drift_snap_last']:+.5f} | {sb['drift_snap_last']:+.5f} |")
add(f"| Drift max abs over the 20 snapshots (%) | {sa['drift_snap_maxabs']:.5f} | {sb['drift_snap_maxabs']:.5f} |")
da, db = sa['drift_mon_vs_first'], sb['drift_mon_vs_first']
add(f"| liquidVolume monitor (gas = {sa['Vtot']:.8e} minus volIntegrate(alpha.water)): drift of last ({da['tL']:.4f} s) vs first ({da['t0']:.4f} s) sample, n = {da['n']} (%) | {da['last']:+.5f} | {db['last']:+.5f} |")
add(f"| Same, max abs over all samples (%) | {da['maxabs']:.5f} | {db['maxabs']:.5f} |")
add(f"| Monitor last vs t = 0 snapshot (%) | {sa['drift_mon_vs_t0']['last']:+.5f} | {sb['drift_mon_vs_t0']['last']:+.5f} |")
add("\nBoth runs gain a trace of gas (monotonic, 0.010% and 0.006% in 4 s); neither is a mass-conservation problem. The two sources agree to 4 digits.\n")
# ---------- per-snapshot
add("## Per-snapshot numbers (both runs)\n")
add("| t (s) | V1a liq | S liq | V1a gas | S gas | V1a jump | S jump (x=2) | S jump (centroid) | V1a xc | S xc | V1a smeared % | S smeared % | V1a mean itf | S mean itf |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")


def fj(x): return 'None' if x is None else f"{x:.3f}"


for ra, rb in zip(a['snaps'], b['snaps']):
    add(f"| {ra['t']:.4f} | {ra['maxliq']:.3f} | {rb['maxliq']:.3f} | {ra['maxgas']:.3f} | {rb['maxgas']:.3f} | {fj(ra['jump'])} | {fj(rb['jump'])} | {fj(rb['jump_c'])} | {ra['xc']:.4f} | {rb['xc']:.4f} | {ra['smeared']:.3f} | {rb['smeared']:.3f} | {ra['meanitf']:.3f} | {rb['meanitf']:.3f} |")
(P3 / 's_rho_result.md').write_text('\n'.join(L) + '\n', encoding='utf-8')
with open(P3 / 'b11_snap.csv', 'w') as f:
    f.write('run,t,max_liquid,max_gas,max_interface,jump_fixed_centre,jump_centroid,Reff,smeared_pct,mean_abs_U_interface,n_interface,gas_volume,centroid_x\n')
    for nm, d in [('V1a_static', a), ('S_rhoST', b)]:
        for r in d['snaps']:
            f.write(f"{nm},{r['t']},{r['maxliq']},{r['maxgas']},{r['maxitf']},{'' if r['jump'] is None else r['jump']},{'' if r['jump_c'] is None else r['jump_c']},{r['Reff']},{r['smeared']},{r['meanitf']},{r['nitf']},{r['gasV']},{r['xc']}\n")
print('written')
