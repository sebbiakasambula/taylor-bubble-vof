"""Check of the grid and time-step findings: Courant margin at T1 and Celik order and GCI for Fr."""
import re, math
import numpy as np
from pathlib import Path
from config import RUNS, WORK
RUN = RUNS
# 1. Courant margin at T1, medium mesh, t >= 5 s
co, aco, dt = [], [], []
t = None
for line in open(RUN / 'V4_medium/log.interFoam', errors='ignore'):
    if line.startswith('Courant Number mean'):
        c = float(line.split('max:')[1])
    elif line.startswith('Interface Courant Number mean'):
        a = float(line.split('max:')[1])
    elif line.startswith('deltaT ='):
        d = float(line.split('=')[1])
    elif line.startswith('Time = '):
        t = float(line.split('=')[1])
        if t >= 5 and 'd' in dir():
            co.append(c); aco.append(a); dt.append(d)
co, aco, dt = map(np.array, (co, aco, dt))
print('steps', len(dt), 'dt median %.4e' % np.median(dt))
print('Co max median %.3f, alphaCo max median %.3f' % (np.median(co), np.median(aco)))
print('ratio Courant-limited dt / actual: by Co 0.5 median %.2f min %.2f; by alphaCo 0.2 median %.2f min %.2f; both median %.2f min %.2f' % (
    np.median(0.5 / co), (0.5 / co).min(), np.median(0.2 / aco), (0.2 / aco).min(),
    np.median(np.minimum(0.5 / co, 0.2 / aco)), np.minimum(0.5 / co, 0.2 / aco).min()))
# 2. Celik order for Fr, 5 s to end
def fit(k, a=5):
    o = np.load(WORK / (k + '.npy')); m = o[:, 0] >= a; return np.polyfit(o[m, 0], o[m, 1], 1)[0]
f1, f2, f3 = fit('V4_fine'), fit('V4_medium'), fit('V4_coarse')
N1, N2, N3 = 90048, 46080, 23667
r21, r32 = math.sqrt(N1 / N2), math.sqrt(N2 / N3)
e21, e32 = f2 - f1, f3 - f2
s = np.sign(e32 / e21)
p = 2.0
for _ in range(200):
    q = math.log((r21 ** p - s) / (r32 ** p - s)); p = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
p_simple = math.log(abs(e32 / e21)) / math.log(r21)
p_14 = math.log(abs(e32 / e21)) / math.log(1.4)
for name, pp, r in (('Celik iterated', p, r21), ('no q term, r21', p_simple, r21), ('r = 1.4', p_14, 1.4)):
    ea = abs(e21 / f1); gf = 1.25 * ea / (r ** pp - 1); gm = 1.25 * ea * r ** pp / (r ** pp - 1)
    print('%-16s p %.4f  GCI fine %.4f%%  GCI medium %.4f%%' % (name, pp, 100 * gf, 100 * gm))
print('Fr %.5f %.5f %.5f  r21 %.4f r32 %.4f' % (f3, f2, f1, r21, r32))
