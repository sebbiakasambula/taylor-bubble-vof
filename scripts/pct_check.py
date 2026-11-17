"""Recompute every comparative percentage in the P2 paper from unrounded data, next to the value the paper states."""
import sys
import numpy as np
from config import WORK
import tb_references as ref
L = lambda k: np.load(WORK / (k + '.npy'))
def fr(t, y, a, b=1e9):
    m = (t >= a) & (t <= b + 1e-3); return np.polyfit(t[m], y[m], 1)[0]
def fm(t, f, a, b=1e9):
    m = (t >= a) & (t <= b + 1e-3); return np.nanmean(f[m])
o = {k: L(k) for k in ['V4_coarse', 'V4_medium', 'V4_fine', 'VA_square', 'T2_prod', 'T3_prod']}
f = {k: L(k) for k in ['f_t1', 'f_t2', 'f_t3', 'f_t1sq']}
OF = lambda k, a, b=1e9: fr(o[k][:, 0], o[k][:, 1], a, b)
FL = lambda k, a, b=1e9: fr(f[k][:, 1], f[k][:, 2], a, b)
v1, v2, v3 = ref.fr_viana(300, 1000), ref.fr_viana(100, 50), ref.fr_viana(20, 500)
l1, l2 = ref.film_llewellin(1000), ref.film_llewellin(50)
pc = lambda a, b: 100 * (a / b - 1)
rows = [
 ('scatter/Fr T1 (15)', 100 * 0.05 / v1), ('scatter/Fr T2 (20)', 100 * 0.05 / v2), ('scatter/Fr T3', 100 * 0.05 / v3),
 ('fine Fr 5-end (0.3267)', OF('V4_fine', 5)), ('fine vs Viana (-3.3)', pc(OF('V4_fine', 5), v1)),
 ('fine film 5-end (0.0659)', fm(o['V4_fine'][:, 0], o['V4_fine'][:, 4], 5)), ('Llewellin T1 (0.0615)', l1),
 ('fine film vs Llewellin (+7.1)', pc(fm(o['V4_fine'][:, 0], o['V4_fine'][:, 4], 5), l1)),
 ('square 5-end (0.3390)', OF('VA_square', 5)), ('square vs fine 5-end (+3.8)', pc(OF('VA_square', 5), OF('V4_fine', 5))),
 ('square 10-end (0.3385)', OF('VA_square', 10)), ('fine 10-end (0.3269)', OF('V4_fine', 10)),
 ('square vs fine 10-end (+3.5)', pc(OF('VA_square', 10), OF('V4_fine', 10))),
 ('square vs medium 5-20 (+3.6)', pc(OF('VA_square', 5, 20), OF('V4_medium', 5, 20))),
 ('Fluent T1 5-20 (0.3376)', FL('f_t1', 5, 20)), ('Fluent T1 vs Viana (-0.1)', pc(FL('f_t1', 5, 20), v1)),
 ('OF T1 5-20 (0.3293)', OF('V4_medium', 5, 20)), ('OF T1 vs Viana', pc(OF('V4_medium', 5, 20), v1)),
 ('Fluent/OF T1 (+2.5)', pc(FL('f_t1', 5, 20), OF('V4_medium', 5, 20))),
 ('Fluent sq vs Fluent str (+0.4)', pc(FL('f_t1sq', 5, 20), FL('f_t1', 5, 20))),
 ('Fluent/OF square (-0.6)', pc(FL('f_t1sq', 5, 20), OF('VA_square', 5, 20))),
 ('OF T2 5-end (0.2180)', OF('T2_prod', 5)), ('OF T2 5-end vs Viana (-12.8)', pc(OF('T2_prod', 5), v2)),
 ('OF T2 film 5-end (0.1345)', fm(o['T2_prod'][:, 0], o['T2_prod'][:, 4], 5)), ('Llewellin T2 (0.139)', l2),
 ('OF T2 film vs Llewellin (-3.3)', pc(fm(o['T2_prod'][:, 0], o['T2_prod'][:, 4], 5), l2)),
 ('OF T3 5-40 (0.3150)', OF('T3_prod', 5, 40)), ('OF T3 5-40 vs Viana (+13.7)', pc(OF('T3_prod', 5, 40), v3)),
 ('Fluent T2 5-20 (0.2637)', FL('f_t2', 5, 20)), ('Fluent T2 vs Viana (+5.5)', pc(FL('f_t2', 5, 20), v2)),
 ('OF T2 5-20 (0.2182)', OF('T2_prod', 5, 20)), ('Fluent/OF T2 (+20.9)', pc(FL('f_t2', 5, 20), OF('T2_prod', 5, 20))),
 ('OF T2 vs Viana 5-20', pc(OF('T2_prod', 5, 20), v2)),
 ('Fluent T2 film 5-20 (0.1422)', fm(f['f_t2'][:, 1], f['f_t2'][:, 4], 5, 20)),
 ('Fluent T2 film vs Llewellin (+2.3)', pc(fm(f['f_t2'][:, 1], f['f_t2'][:, 4], 5, 20), l2)),
 ('OF T2 film 5-20', fm(o['T2_prod'][:, 0], o['T2_prod'][:, 4], 5, 20)),
 ('Fluent/OF T2 film (+5.7)', pc(fm(f['f_t2'][:, 1], f['f_t2'][:, 4], 5, 20), fm(o['T2_prod'][:, 0], o['T2_prod'][:, 4], 5, 20))),
 ('Fluent T3 5-20 (0.2829)', FL('f_t3', 5, 20)), ('Fluent T3 vs Viana (+2.1)', pc(FL('f_t3', 5, 20), v3)),
 ('OF T3 5-20 (0.3170)', OF('T3_prod', 5, 20)), ('Fluent/OF T3 (-10.8)', pc(FL('f_t3', 5, 20), OF('T3_prod', 5, 20))),
 ('OF T3 vs Viana 5-20 (+14.4)', pc(OF('T3_prod', 5, 20), v3)),
 ('T1 cross-code / GCImed 2.33 (1.1)', abs(pc(FL('f_t1', 5, 20), OF('V4_medium', 5, 20))) / 2.33),
 ('T2 cross-code / 2.33', abs(pc(FL('f_t2', 5, 20), OF('T2_prod', 5, 20))) / 2.33),
 ('T3 cross-code / 2.33', abs(pc(FL('f_t3', 5, 20), OF('T3_prod', 5, 20))) / 2.33),
]
for k, v in rows:
    print(f"{k:38s} {v:9.4f}")
