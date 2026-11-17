"""Set 17, part 3. The same film check as set17_profiles.py at T1 (Nf 1000, window from 2.0 D below the nose)."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import set17_profiles as P  # noqa: E402
from tb_params import fluid  # noqa: E402

f1 = fluid(300, 1000, 1000)
P.G = (f1["rho_l"] - f1["rho_g"]) / f1["mu_l"]
P.CLEAR = 2.0


def folders(name, lo, hi):
    out = []
    for d in (P.RUNS / name).iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if lo <= t <= hi:
            out.append((t, d.name))
    return [n for _, n in sorted(out)]


lines = ["# Set 17, part 3. Film velocity at T1 against Brown's profile", "",
         "G = %.4f, window from 2.0 D below the nose to 0.5 D above the axis tail, snapshots in 5 to 10 s." % P.G, "",
         "| run | snapshots | columns | film flux / theory, median (10th to 90th) | interface-cell speed / theory, median (10th to 90th) | mean film in these columns | columns by distance behind the nose (D) |",
         "|---|---|---|---|---|---|---|"]
for lab, name in (("OF baseline", "V4_medium"), ("OF predictor", "T1med_mompred"), ("OF geometric and predictor", "T1med_A2plic")):
    ts = folders(name, 5.0, 10.0)
    q_all, i_all, d_all, xs = [], [], [], []
    for t in ts:
        x, y, al, ux, wall = P.of_snap(name, t)
        nose, cols = P.columns(x, y, al, ux, wall)
        q, i = P.stats(cols, wall)
        q_all += list(q); i_all += list(i); d_all += [c[1] for c in cols]; xs += [-c[0] for c in cols]
    q_all, i_all = np.array(q_all), np.array(i_all)
    lines.append("| %s | %d | %d | %.4f (%.4f to %.4f) | %.4f (%.4f to %.4f) | %.4f | %.2f to %.2f |" % (
        lab, len(ts), len(q_all), np.median(q_all), *np.percentile(q_all, [10, 90]), np.median(i_all),
        *np.percentile(i_all, [10, 90]), np.mean(d_all), min(xs), max(xs)))
steps = [5850, 6435, 7020, 7605, 8190, 8775, 9360, 9945]
q_all, i_all, d_all, xs = [], [], [], []
for s in steps:
    x, y, al, ux, wall = P.fl_snap("T1_Case", "s_snap", s)
    nose, cols = P.columns(x, y, al, ux, wall)
    q, i = P.stats(cols, wall)
    q_all += list(q); i_all += list(i); d_all += [c[1] for c in cols]; xs += [-c[0] for c in cols]
q_all, i_all = np.array(q_all), np.array(i_all)
lines.append("| Fluent | %d | %d | %.4f (%.4f to %.4f) | %.4f (%.4f to %.4f) | %.4f | %.2f to %.2f |" % (
    len(steps), len(q_all), np.median(q_all), *np.percentile(q_all, [10, 90]), np.median(i_all),
    *np.percentile(i_all, [10, 90]), np.mean(d_all), min(xs), max(xs)))
(HERE / "set17_profiles_t1_out.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
