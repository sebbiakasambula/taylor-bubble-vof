"""B11: V1a_static (standard interFoam) against S_rhoST (density-weighted surface tension), same mesh.
Uses scal/vec and tb_post exactly as snap.py does (read from snap.py itself, as tovtk.py does)."""
import re, sys, json
from pathlib import Path
import numpy as np
from config import RUNS, WORK
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])   # gives scal, vec, tb_post, np, re, Path
RUN = RUNS
SIG = 0.0999
OUT = {}

def monitor(case):
    t = []; m = []; xl = []; yl = []
    for l in open(RUN / case / 'postProcessing/maxSpeed/0/fieldMinMax.dat'):
        if l.startswith('#') or not l.strip(): continue
        g = re.match(r'\s*([0-9.e+-]+)\s+mag\(U\)\s+([0-9.e+-]+)\s+\([^)]*\)\s+([0-9.e+-]+)\s+\(\s*([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)\s*\)', l)
        t.append(float(g[1])); m.append(float(g[3])); xl.append(float(g[4])); yl.append(float(g[5]))
    return np.array(t), np.array(m), np.array(xl), np.array(yl)

def tdirs(c):
    return sorted([d for d in c.iterdir() if d.is_dir() and re.fullmatch(r'[0-9.]+', d.name) and d.name != '0'], key=lambda d: float(d.name))

def smeared(a, V):
    g = (1 - a) * V
    return 100 * g[a > 0.5].sum() / g.sum()

def run(case):
    c = RUN / case
    C = vec(c / '0/C'); V = scal(c / '0/V'); x, y = C[:, 0], C[:, 1]
    a0 = scal(c / '0/alpha.water')
    R = {'case': case, 'ncell': len(V)}
    # --- t = 0 shape
    R['t0'] = dict(Reff=tb_post.effective_radius(a0, V), gasV=tb_post.gas_volume(a0, V), smeared=smeared(a0, V))
    # --- snapshots
    rows = []
    for d in tdirs(c):
        t = float(d.name)
        a = scal(d / 'alpha.water'); U = vec(d / 'U'); p = scal(d / 'p'); sp = np.linalg.norm(U, axis=1)
        gas = a < 0.01; liq = a > 0.99; itf = ~gas & ~liq
        Reff = tb_post.effective_radius(a, V)
        j = tb_post.laplace_jump(x, y, a, p, 2.0, 0.25)
        xc = tb_post.centroid_x(x, a, V)
        jc = tb_post.laplace_jump(x, y, a, p, xc, 0.25)
        gm = (1 - a) * V; m5 = a < 0.5
        rows.append(dict(t=t, xc=xc, jump_c=jc, xmin5=float(x[m5].min()), xmax5=float(x[m5].max()), rmax5=float(y[m5].max()),
                         n_in=int(((np.hypot(x - 2.0, y) < 0.175) & gas).sum()), maxgas=sp[gas].max(), maxitf=sp[itf].max(), maxliq=sp[liq].max(), maxall=sp.max(),
                         gasV=tb_post.gas_volume(a, V), Reff=Reff, jump=j, target=2 * SIG / Reff,
                         smeared=smeared(a, V), meanitf=sp[itf].mean(), nitf=int(itf.sum()),
                         meanitf_vw=float(np.sum(sp[itf] * V[itf]) / V[itf].sum())))
    R['snaps'] = rows
    # --- last snapshot: where is the maximum
    d = tdirs(c)[-1]
    a = scal(d / 'alpha.water'); U = vec(d / 'U'); sp = np.linalg.norm(U, axis=1)
    k = int(np.argmax(sp)); al = a[k]
    cls = 'gas (alpha<0.01)' if al < 0.01 else ('liquid (alpha>0.99)' if al > 0.99 else 'interface (0.01<=alpha<=0.99)')
    R['last'] = dict(t=float(d.name), max=float(sp[k]), x=float(x[k]), r=float(y[k]), alpha=float(al), cls=cls,
                     Ux=float(U[k, 0]), Ur=float(U[k, 1]))
    # --- monitor
    t, m, xl, yl = monitor(case)
    R['mon'] = dict(t=t, m=m, xl=xl, yl=yl)
    # --- liquid volume monitor
    hdr = open(c / 'postProcessing/liquidVolume/0/volFieldValue.dat').read().splitlines()
    Vtot = float([l for l in hdr if l.startswith('# Volume')][0].split(':')[1])
    dd = np.array([[float(q) for q in l.split()] for l in hdr if not l.startswith('#') and l.strip()])
    R['lv'] = dict(t=dd[:, 0], liq=dd[:, 1], Vtot=Vtot)
    return R

def rng(v): return (float(np.nanmin(v)), float(np.nanmax(v)))

def summarise(R):
    s = {}
    t, m = R['mon']['t'], R['mon']['m']
    k = (t >= 0.3) & (t <= 1.0001)
    s['ms_1s'] = float(m[k].mean()); s['ms_1s_n'] = int(k.sum()); s['ms_1s_max'] = float(m[k].max())
    s['ms_1s_trap'] = float(np.trapezoid(m[k], t[k]) / (t[k][-1] - t[k][0]))
    k = (t >= 0.3)
    s['ms_all'] = float(m[k].mean()); s['ms_all_n'] = int(k.sum()); s['ms_all_max'] = float(m[k].max())
    s['ms_all_trap'] = float(np.trapezoid(m[k], t[k]) / (t[k][-1] - t[k][0]))
    s['ms_peak_run'] = float(m.max()); s['t_peak_run'] = float(t[np.argmax(m)])
    s['mon_t_first'] = float(t[k][0]); s['mon_t_last'] = float(t[-1]); s['mon_last'] = float(m[-1])
    s['mon_last_loc'] = (float(R['mon']['xl'][-1]), float(R['mon']['yl'][-1]))
    # windowed means of the monitor
    edges = [0.3, 1.0, 2.0, 3.0, 3.96]
    s['mon_win'] = [(edges[i], edges[i + 1], float(m[(t >= edges[i]) & (t < edges[i + 1])].mean())) for i in range(len(edges) - 1)]
    sn = R['snaps']
    T = np.array([r['t'] for r in sn])
    def col(n): return np.array([np.nan if r[n] is None else r[n] for r in sn])
    for name, sel in [('1s', T <= 1.0 + 1e-3), ('all20', T > 0), ('ge03', T >= 0.3)]:
        s['n_' + name] = int(sel.sum())
        for q in ['maxliq', 'maxgas', 'maxitf', 'jump', 'target']:
            v = col(q)[sel]
            s[f'{q}_{name}'] = dict(min=float(np.nanmin(v)), max=float(np.nanmax(v)), mean=float(np.nanmean(v)), sd=float(np.nanstd(v, ddof=1)) if len(v) > 1 else None)
        jj = col('jump')[sel]; tt = col('target')[sel]; ok = ~np.isnan(jj)
        s['jump_nvalid_' + name] = int(ok.sum())
        s['jump_err_' + name] = float(100 * (jj[ok].mean() / tt[ok].mean() - 1))
        jc = np.array([np.nan if r['jump_c'] is None else r['jump_c'] for r in sn])[sel]; okc = ~np.isnan(jc)
        s['jumpc_' + name] = dict(n=int(okc.sum()), min=float(jc[okc].min()), max=float(jc[okc].max()), mean=float(jc[okc].mean()), sd=float(jc[okc].std(ddof=1)), err=float(100 * (jc[okc].mean() / tt[okc].mean() - 1)))
    # shape
    s['shape'] = dict(Reff0=R['t0']['Reff'], Reff_first=sn[0]['Reff'], Reff_last=sn[-1]['Reff'],
                      sm0=R['t0']['smeared'], sm_first=sn[0]['smeared'], sm_last=sn[-1]['smeared'],
                      Reff_min=float(min(r['Reff'] for r in sn)), Reff_max=float(max(r['Reff'] for r in sn)))
    s['xc_series'] = [(r['t'], r['xc']) for r in sn]
    s['xc_last'] = sn[-1]['xc']; s['dx_last'] = sn[-1]['xc'] - 2.0
    tt_ = np.array([r['t'] for r in sn]); xx_ = np.array([r['xc'] for r in sn])
    kk = tt_ >= 2.3
    s['xc_slope_ge23'] = float(np.polyfit(tt_[kk], xx_[kk], 1)[0])
    s['xc_slope_all'] = float(np.polyfit(tt_, xx_, 1)[0])
    s['xc_range'] = (float(xx_.min()), float(xx_.max()))
    s['ext_x_range5'] = (min(r['xmin5'] for r in sn), max(r['xmax5'] for r in sn))
    s['ext_r_range5'] = (min(r['rmax5'] for r in sn), max(r['rmax5'] for r in sn))
    s['meanitf_all'] = float(np.mean([r['meanitf'] for r in sn])); s['nitf_range'] = (min(r['nitf'] for r in sn), max(r['nitf'] for r in sn))
    s['nin_min'] = min(r['n_in'] for r in sn)
    s['last'] = R['last']
    s['last_snap'] = dict(maxgas=sn[-1]['maxgas'], maxitf=sn[-1]['maxitf'], maxliq=sn[-1]['maxliq'], meanitf=sn[-1]['meanitf'], nitf=sn[-1]['nitf'], meanitf_vw=sn[-1]['meanitf_vw'])
    s['meanitf_t'] = [(r['t'], r['meanitf']) for r in sn]
    # drift
    gas0_snap = R['t0']['gasV']
    gs = np.array([r['gasV'] for r in sn])
    s['drift_snap_last'] = float(100 * (gs[-1] / gas0_snap - 1)); s['drift_snap_maxabs'] = float(100 * np.abs(gs / gas0_snap - 1).max())
    gas_mon = R['lv']['Vtot'] - R['lv']['liq']
    s['drift_mon_vs_first'] = dict(n=len(gas_mon), gas0=float(gas_mon[0]), gasL=float(gas_mon[-1]), last=float(100 * (gas_mon[-1] / gas_mon[0] - 1)), maxabs=float(100 * np.abs(gas_mon / gas_mon[0] - 1).max()), t0=float(R['lv']['t'][0]), tL=float(R['lv']['t'][-1]))
    s['drift_mon_vs_t0'] = dict(last=float(100 * (gas_mon[-1] / gas0_snap - 1)), maxabs=float(100 * np.abs(gas_mon / gas0_snap - 1).max()), first=float(100 * (gas_mon[0] / gas0_snap - 1)))
    s['gas0_snap'] = gas0_snap; s['gasV_last'] = float(gs[-1]); s['Vtot'] = R['lv']['Vtot']
    return s

if __name__ == '__main__':
    res = {}
    for case in ['V1a_static', 'S_rhoST']:
        R = run(case); res[case] = (R, summarise(R))
    import pickle
    pickle.dump({k: (dict(snaps=v[0]['snaps'], t0=v[0]['t0'], last=v[0]['last']), v[1]) for k, v in res.items()}, open(WORK / 'b11_res.pkl', 'wb'))
    for case, (R, s) in res.items():
        print('=====', case, 'cells', R['ncell'])
        print('per snapshot')
        for r in R['snaps']:
            print(f"xc {r['xc']:.4f} jc {r['jump_c'] if r['jump_c'] is None else round(r['jump_c'],4)} t {r['t']:7.4f} liq {r['maxliq']:.4f} gas {r['maxgas']:.3f} itf {r['maxitf']:.3f} jump {r['jump'] if r['jump'] is None else round(r['jump'],4)} Reff {r['Reff']:.5f} smeared {r['smeared']:.3f} meanitf {r['meanitf']:.4f} gasV {r['gasV']:.8e}")
        for k, v in s.items():
            if k not in ('meanitf_t',): print(' ', k, v)
