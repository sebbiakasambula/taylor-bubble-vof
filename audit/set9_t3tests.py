"""Data audit, set 9: tests at T3. OpenFOAM T3_A5_plic (A2 settings: interIsoFoam plicRDF + momentum predictor), read
with the body nose (of2/a5body.py, rerun into audit/a5) because its axis cells are unreliable; Fluent D1 (Green-Gauss
node based gradient) and D2 (Modified HRIC) from the per-snapshot tables built from the raw exports (fl2/run_new.py,
rebuilt in audit/flcopy). Window 5 to 8.4 s. Not in the paper yet. Writes set9_t3tests_out.md.
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FR_OF, FR_F = 0.320527, 0.282005
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3


def f4(v):
    return "%.4f" % v


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    n = len(t)
    return float(p[0]), (math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)) if n > 2 else float("nan")), n


def load(p):
    rows = list(csv.DictReader(open(p)))
    M = {}
    for k in rows[0]:
        try:
            M[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            M[k] = np.array([r[k] for r in rows])
    return M


out = ["## Set 9. Tests at T3 (not in the paper yet), recomputed 14 Nov\n"]
A5 = np.load(HERE / "a5" / "T3_A5_plic.npy")
P = np.load(HERE / "a5" / "T3_prod.npy")
fa = fit(A5[:, 0], A5[:, 1])
idx = [int(np.argmin(np.abs(A5[:, 0] - x))) for x in P[:, 0]]
fm = fit(A5[idx, 0], A5[idx, 1])
M5, MP = load(HERE / "m" / "T3_A5_plic.csv"), load(HERE / "m" / "T3_prod.csv")
w5 = (M5["t"] >= LO) & (M5["t"] <= HI)
wp = (MP["t"] >= LO) & (MP["t"] <= HI)
out += ["### OpenFOAM A5 (A2 settings at T3), body nose\n",
        "- Fr all %s (s.e. %s, n %d), %s%%, phi %s; matched (nearest the 5 T3_prod times) %s (n %d), %s%%, phi %s" % (
            f4(fa[0]), "%.6f" % fa[1], fa[2], f4(100 * (fa[0] / FR_OF - 1)), f4((fa[0] - FR_OF) / (FR_F - FR_OF)),
            f4(fm[0]), fm[2], f4(100 * (fm[0] / FR_OF - 1)), f4((fm[0] - FR_OF) / (FR_F - FR_OF))),
        "- T3_prod body nose %s (identical to its axis nose 0.320527); baseline s.e. 0.010495 (3.27%% of Fr)" % f4(fit(P[:, 0], P[:, 1])[0]),
        "- Fluent against A5: all %s%%, matched %s%%; against the baseline %s%%" % (
            f4(100 * (FR_F / fa[0] - 1)), f4(100 * (FR_F / fm[0] - 1)), f4(100 * (FR_F / FR_OF - 1))),
        "- gas outside the main body: A5 mean %s%%, T3_prod %s%%" % (f4(A5[:, 3].mean()), f4(P[:, 3].mean())),
        "- smeared gas (metrics smear, independent of the nose): A5 %s%%, T3_prod %s%%" % (f4(M5["smear"][w5].mean()), f4(MP["smear"][wp].mean())),
        "- median largest gas speed: A5 %s, T3_prod %s m/s" % (f4(np.median(M5["pk_gas"][w5])), f4(np.median(MP["pk_gas"][wp]))),
        "- A5 largest gas drift %s%%" % f4(100 * np.max(np.abs(M5["G"] / M5["G"][0] - 1)))]

FB = load(HERE / "m" / "F_T3.csv")
wb = (FB["t"] >= LO) & (FB["t"] <= HI)
tb = FB["t"][wb]
base = fit(FB["t"][wb], FB["nose"][wb])[0]
out += ["\n### Fluent D1 and D2\n", "Fluent T3 baseline over 5 to 8.4 s: Fr %s (n %d), film %s, out_b %s%%, smear %s%%.\n" % (
    f4(base), int(wb.sum()), f4(FB["film"][wb].mean()), f4(FB["out_b"][wb].mean()), f4(FB["smear"][wb].mean())),
        "| run | Fr all (n) | Fr matched | shift all (%) | shift matched (%) | phi_F all | paired slope of the difference (%) | film | out_b (%) | smear (%) |",
        "|---|---|---|---|---|---|---|---|---|---|"]
for n, lab in (("T3D1", "D1 Green-Gauss node based gradient"), ("T3D2", "D2 Modified HRIC")):
    p = HERE / "flcopy" / ("m_F_T3_%s.csv" % n[2:])
    if not p.exists():
        p = HERE / "flcopy" / ("m_F_T3_%s.csv" % n[2:])
    M = load(p)
    w = (M["t"] >= LO) & (M["t"] <= HI) & np.isfinite(M["nose"])
    a = fit(M["t"][w], M["nose"][w])
    idx = [int(np.argmin(np.abs(M["t"] - x))) for x in tb]
    m = fit(M["t"][idx], M["nose"][idx])[0]
    d = fit(tb, M["nose"][idx] - FB["nose"][wb])
    out.append("| %s | %s (%d) | %s | %s | %s | %s | %s +- %s | %s | %s | %s |" % (
        lab, f4(a[0]), a[2], f4(m), f4(100 * (a[0] / base - 1)), f4(100 * (m / base - 1)), f4((base - a[0]) / (base - FR_OF)),
        f4(100 * d[0] / base), f4(100 * d[1] / base), f4(M["film"][w].mean()), f4(M["out_b"][w].mean()), f4(M["smear"][w].mean())))
    for hi in (7.79, 7.11):
        wh = (M["t"] >= LO) & (M["t"] <= hi) & np.isfinite(M["nose"])
        wbh = (FB["t"] >= LO) & (FB["t"] <= hi)
        bh = fit(FB["t"][wbh], FB["nose"][wbh])[0]
        ih = [int(np.argmin(np.abs(M["t"] - x))) for x in FB["t"][wbh]]
        out.append("  window to %s s: all %s%%, matched %s%%" % (hi, f4(100 * (fit(M["t"][wh], M["nose"][wh])[0] / bh - 1)),
                                                           f4(100 * (fit(M["t"][ih], M["nose"][ih])[0] / bh - 1))))
(HERE / "set9_t3tests_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
