"""Data audit, set 15: T2_M1c_mompred_coarse, the momentum predictor (MULES, 1 outer loop) on the coarse T2 mesh (hint 2).
Compared with the coarse baseline T2_C1_coarse, the medium predictor run T2_mompred, the medium baseline T2_prod, and A2 on
the coarse and medium meshes. Window 5 to 8.4 s. 'all' = every snapshot in the window; 'matched' = the snapshot nearest
each of the four T2_prod window times (as sets 5 and 7). Fluent 0.263163 (medium). Writes set15_m1c_out.md.
"""
import csv
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
LO, HI = 5.0, 8.4
FR_F = 0.263163


def f4(v):
    return "%.4f" % v


def load(name):
    rows = list(csv.DictReader(open(HERE / "m" / ("%s.csv" % name))))
    M = {}
    for k in rows[0]:
        try:
            M[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            M[k] = np.array([r[k] for r in rows])
    return M


def win(M):
    return (M["t"] >= LO - 1e-6) & (M["t"] <= HI + 1e-3)


def fit(t, x):
    if len(t) > 3:
        p, cov = np.polyfit(t, x, 1, cov=True)
        return float(p[0]), float(np.sqrt(cov[0, 0])), len(t)
    return float(np.polyfit(t, x, 1)[0]), float("nan"), len(t)


BASE = load("T2_prod")
TB = BASE["t"][win(BASE) & np.isfinite(BASE["nose"])]


def fr_all(M):
    w = win(M) & np.isfinite(M["nose"])
    return fit(M["t"][w], M["nose"][w])


def fr_matched(M):
    idx = [int(np.argmin(np.abs(M["t"] - t))) for t in TB]
    return fit(M["t"][idx], M["nose"][idx])[0], float(np.max(np.abs(M["t"][idx] - TB)))


def mean(M, k):
    w = win(M) & np.isfinite(M[k])
    return float(M[k][w].mean())


def median(M, k):
    w = win(M) & np.isfinite(M[k])
    return float(np.median(M[k][w]))


def wall(name):
    ct = re.findall(r"ClockTime = ([0-9.]+) s", (RUNS / name / "log.interFoam").read_text(errors="replace"))
    return float(ct[-1]) / 3600


def pct(a, b):
    return 100 * (a / b - 1)


names = [("T2_C1_coarse", "baseline, coarse"), ("T2_prod", "baseline, medium"), ("T2_M1c_mompred_coarse", "predictor, coarse"),
         ("T2_mompred", "predictor, medium"), ("T2_A2c_coarse", "A2, coarse"), ("T2_A2_plic", "A2, medium")]
R = {}
out = ["## Set 15. T2_M1c_mompred_coarse, the momentum predictor on the coarse T2 mesh (hint 2)\n",
       "| run | what | Fr all (s.e., n) | Fr matched (max offset s) | film | smear (%) | L_b | out_b (%) | pk_all median | pk_liq median |",
       "|---|---|---|---|---|---|---|---|---|---|"]
for nm, lab in names:
    M = load(nm)
    a, (m, off) = fr_all(M), fr_matched(M)
    R[nm] = dict(a=a[0], m=m, film=mean(M, "film"), M=M)
    out.append("| %s | %s | %s (%s, %d) | %s (%.4f) | %s | %s | %s | %s | %s | %s |" % (
        nm, lab, f4(a[0]), f4(a[1]), a[2], f4(m), off, f4(mean(M, "film")), f4(mean(M, "smear")), f4(mean(M, "L_b")),
        f4(mean(M, "out_b")), f4(median(M, "pk_all")), f4(median(M, "pk_liq"))))
out.append("")
for (c, mm, lab) in (("T2_C1_coarse", "T2_prod", "baseline"), ("T2_M1c_mompred_coarse", "T2_mompred", "predictor"),
                     ("T2_A2c_coarse", "T2_A2_plic", "A2")):
    out.append("- %s, medium against coarse: Fr all %s%%, matched %s%%; film all %s%%" % (
        lab, f4(pct(R[mm]["a"], R[c]["a"])), f4(pct(R[mm]["m"], R[c]["m"])), f4(pct(R[mm]["film"], R[c]["film"]))))
for mesh, b, p in (("coarse", "T2_C1_coarse", "T2_M1c_mompred_coarse"), ("medium", "T2_prod", "T2_mompred")):
    out.append("- predictor against baseline on the %s mesh: Fr all %s%%, matched %s%%; film %s%%" % (
        mesh, f4(pct(R[p]["a"], R[b]["a"])), f4(pct(R[p]["m"], R[b]["m"])), f4(pct(R[p]["film"], R[b]["film"]))))
M = R["T2_M1c_mompred_coarse"]["M"]
g = M["G"]
out += ["- M1c gas drift, last against t = 0: %s%%; starting gas against T2_C1_coarse: %s%%" % (
            f4(pct(g[-1], g[0])), f4(pct(g[0], R["T2_C1_coarse"]["M"]["G"][0]))),
        "- M1c end time %s s, %d snapshots, wall time %s h" % (f4(M["t"][-1]), len(M["t"]), f4(wall("T2_M1c_mompred_coarse")))]
(HERE / "set15_m1c_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
