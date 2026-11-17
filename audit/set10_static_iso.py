"""Data audit, set 10: static isoAdvector, curvature-smoothing and normal-gradient tests on case S (not in the paper).
Runs: S_B2_plic, S_B2b_plicpred, S_B2c_isopred, S_B2d_plic10, S_B2e_mules10, S_B3a_mules_sm2, S_B3b_plic_sm2,
S_G1a_mules_pl, S_G1b_plic_pl, with V1a_static as the MULES reference. Definitions as set 1 and b2chk.py: mean maximum
speed = sample mean of the maxSpeed monitor over 0.3 to 1 s; liquid peak, gas peak, jump and Laplace value from the
written snapshots with 0 < t <= 1.0001 s; wedge half-angle read from the mesh. Writes set10_static_iso_out.md.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import tb_post  # noqa: E402

sys.path.insert(0, str(HERE))
from of_metrics import scal, vec, half_angle, RUNS  # noqa: E402


def f4(v):
    return "nan" if v is None or (isinstance(v, float) and math.isnan(v)) else "%.4f" % v


def sigma(case):
    return float(re.search(r"sigma\s+([0-9.eE+-]+)", (case / "constant/transportProperties").read_text())[1])


def app(case):
    return re.search(r"application\s+(\w+)", (case / "system/controlDict").read_text())[1]


out = ["## Set 10. Static isoAdvector and related tests on case S (not in the paper), recomputed 14 Nov\n",
       "| run | solver | wedge half-angle | monitor samples 0.3-1 s | mean max speed (m/s) | max in 0.3-1 s | snapshots | liquid peak min | liquid peak max | gas peak max | jump min | jump max | jump mean | Laplace | error (%) |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for name in ("V1a_static", "S_B2_plic", "S_B2b_plicpred", "S_B2c_isopred", "S_B2d_plic10", "S_B2e_mules10", "S_B3a_mules_sm2",
             "S_B3b_plic_sm2", "S_G1a_mules_pl", "S_G1b_plic_pl"):
    case = RUNS / name
    ha = half_angle(case)
    t, m = [], []
    for line in open(case / "postProcessing/maxSpeed/0/fieldMinMax.dat"):
        g = re.match(r"\s*([0-9.e+-]+)\s+mag\(U\)\s+([0-9.e+-]+)\s+\([^)]*\)\s+([0-9.e+-]+)", line)
        if g:
            t.append(float(g[1]))
            m.append(float(g[3]))
    t, m = np.array(t), np.array(m)
    k = (t >= 0.3) & (t <= 1.0001)
    C = vec(case / "0/C")
    V = scal(case / "0/V")
    sg = sigma(case)
    liq, gas, jp, lap = [], [], [], []
    for d in sorted([d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name) and d.name != "0"], key=lambda d: float(d.name)):
        if float(d.name) > 1.0001:
            continue
        a = scal(d / "alpha.water")
        U = vec(d / "U")
        p = scal(d / "p")
        sp = np.linalg.norm(U, axis=1)
        liq.append(sp[a > 0.99].max())
        gas.append(sp[a < 0.01].max())
        j = tb_post.laplace_jump(C[:, 0], C[:, 1], a, p, 2.0, 0.25)
        jp.append(np.nan if j is None else j)
        lap.append(2 * sg / tb_post.effective_radius(a, V, ha))
    jp = np.array(jp)
    out.append("| %s | %s | %s | %d | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        name, app(case), f4(ha), int(k.sum()), f4(float(m[k].mean())), f4(float(m[k].max())), len(liq), f4(min(liq)), f4(max(liq)),
        f4(max(gas)), f4(float(np.nanmin(jp))), f4(float(np.nanmax(jp))), f4(float(np.nanmean(jp))), f4(float(np.mean(lap))),
        f4(100 * (float(np.nanmean(jp)) / float(np.mean(lap)) - 1))))
    print("done", name, file=sys.stderr)
(HERE / "set10_static_iso_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
