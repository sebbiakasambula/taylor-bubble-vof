"""Data audit, set 6: momentum predictor (T2_mompred, T1med_mompred, T1sq_mompred) and density-weighted surface tension
(T2_rhoST, T3_rhoST), plus the build check of the rebuilt solver (T2_rhoST_check against T2_prod, first 64 steps).
Per-snapshot tables from audit/m (recomputed). Writes set6_tests_out.md.
"""
import csv
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
FR_F_T2 = 0.263163
VIANA_T2 = 0.250410


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


def win(M, lo, hi):
    return (M["t"] >= lo - 1e-6) & (M["t"] <= hi + 1e-3)


def fr(M, lo, hi):
    w = win(M, lo, hi) & np.isfinite(M["nose"])
    p = np.polyfit(M["t"][w], M["nose"][w], 1)
    return float(p[0]), int(w.sum())


def mean(M, k, lo, hi):
    w = win(M, lo, hi) & np.isfinite(M[k])
    return float(M[k][w].mean()), int(w.sum())


def wall(name):
    ct = re.findall(r"ClockTime = ([0-9.]+) s", (RUNS / name / "log.interFoam").read_text(errors="replace"))
    return float(ct[-1]) / 3600


out = ["## Set 6. Momentum predictor and density weighting, recomputed 14 Nov\n"]
T2, MP = load("T2_prod"), load("T2_mompred")
b, m = fr(T2, 5, 8.4), fr(MP, 5, 8.4)
gap_b, gap_m = 100 * (FR_F_T2 / b[0] - 1), 100 * (FR_F_T2 / m[0] - 1)
out += ["### Momentum predictor at T2 (5 to 8.4 s)\n",
        "- Fr without %s (n %d), with %s (n %d): %s%%; gap to Fluent (0.263163) %s%% to %s%%; phi %s (closed share of the gap %s)" % (
            f4(b[0]), b[1], f4(m[0]), m[1], f4(100 * (m[0] / b[0] - 1)), f4(gap_b), f4(gap_m), f4((m[0] - b[0]) / (FR_F_T2 - b[0])),
            f4(1 - gap_m / gap_b)),
        "- film without %s, with %s" % (f4(mean(T2, "film", 5, 8.4)[0]), f4(mean(MP, "film", 5, 8.4)[0])),
        "- wall time T2_mompred %s h" % f4(wall("T2_mompred"))]

V4, T1m, VA, T1s = load("V4_medium"), load("T1med_mompred"), load("VA_square"), load("T1sq_mompred")
F1 = load("F_T1")
fF = fr(F1, 5, 10)
a, c = fr(V4, 5, 10), fr(T1m, 5, 10)
s10, s10m = fr(VA, 10, 15), fr(T1s, 10, 15)
s5, s5m = fr(VA, 5, 10), fr(T1s, 5, 10)
out += ["\n### Momentum predictor at T1\n",
        "- medium 5 to 10 s: without %s, with %s, %s%%; Fluent 5 to 10 s %s; gap %s%% to %s%%" % (
            f4(a[0]), f4(c[0]), f4(100 * (c[0] / a[0] - 1)), f4(fF[0]), f4(100 * (fF[0] / a[0] - 1)), f4(100 * (fF[0] / c[0] - 1))),
        "- medium film 5 to 10 s: without %s, with %s" % (f4(mean(V4, "film", 5, 10)[0]), f4(mean(T1m, "film", 5, 10)[0])),
        "- square 10 to 15 s: without %s, with %s, %s%%" % (f4(s10[0]), f4(s10m[0]), f4(100 * (s10m[0] / s10[0] - 1))),
        "- square 5 to 10 s: without %s, with %s, %s%%" % (f4(s5[0]), f4(s5m[0]), f4(100 * (s5m[0] / s5[0] - 1))),
        "- square against medium over 5 to 10 s: with the predictor %s%%, without %s%%" % (
            f4(100 * (s5m[0] / c[0] - 1)), f4(100 * (s5[0] / a[0] - 1))),
        "- wall time T1med_mompred %s h, T1sq_mompred %s h; end times %s s and %s s" % (
            f4(wall("T1med_mompred")), f4(wall("T1sq_mompred")), f4(T1m["t"][-1]), f4(T1s["t"][-1]))]

RH = load("T2_rhoST")
r = fr(RH, 5, 8.4)
fb, fr_ = mean(T2, "film", 5, 8.4)[0], mean(RH, "film", 5, 8.4)[0]
lb, lr = mean(T2, "L_b", 5, 8.4)[0], mean(RH, "L_b", 5, 8.4)[0]
out += ["\n### Density weighting at T2 (5 to 8.4 s)\n",
        "- Fr without %s, with %s, %s%%; gap to Fluent %s%% to %s%%; against Viana %s%%" % (
            f4(b[0]), f4(r[0]), f4(100 * (r[0] / b[0] - 1)), f4(gap_b), f4(100 * (FR_F_T2 / r[0] - 1)), f4(100 * (r[0] / VIANA_T2 - 1))),
        "- film %s to %s (%s%%); bubble length %s to %s (%s%%)" % (f4(fb), f4(fr_), f4(100 * (fr_ / fb - 1)), f4(lb), f4(lr), f4(100 * (lr / lb - 1))),
        "- largest gas drift %s%%; out_b mean %s%%" % (f4(100 * np.max(np.abs(RH["G"] / RH["G"][0] - 1))), f4(mean(RH, "out_b", 5, 8.4)[0]))]

RT = load("T3_rhoST")
T3 = load("T3_prod")
w3, w3r = win(T3, 5, 8.4), win(RT, 5, 8.4)
out += ["\n### Density weighting at T3, quantities from the per-snapshot tables\n",
        "- smeared gas over 5 to 8.4 s: without %s to %s%%, with %s to %s%%" % (
            f4(T3["smear"][w3].min()), f4(T3["smear"][w3].max()), f4(RT["smear"][w3r].min()), f4(RT["smear"][w3r].max())),
        "- median largest speed 5 to 8.4 s: without %s, with %s m/s" % (f4(np.median(T3["pk_all"][w3])), f4(np.median(RT["pk_all"][w3r]))),
        "- largest gas drift with the weight %s%%" % f4(100 * np.max(np.abs(RT["G"] / RT["G"][0] - 1)))]
for tt in (1.0, 2.2, 2.4, 5.0, 8.4):
    k = int(np.argmin(np.abs(RT["t"] - tt)))
    out.append("- with the weight at t = %s: largest speed %s, largest liquid speed %s, nose %s" % (f4(RT["t"][k]), f4(RT["pk_all"][k]), f4(RT["pk_liq"][k]), f4(RT["nose"][k])))
for tt in (2.26, 5.0):
    k = int(np.argmin(np.abs(T3["t"] - tt)))
    out.append("- without the weight at t = %s: largest speed %s, largest liquid speed %s, nose %s" % (f4(T3["t"][k]), f4(T3["pk_all"][k]), f4(T3["pk_liq"][k]), f4(T3["nose"][k])))
out.append("- nose at t = 0: without %s, with %s" % (f4(T3["nose"][0]), f4(RT["nose"][0])))

# build check: T2_rhoST_check against T2_prod, step by step
def steps(name, n):
    rows, cur = [], {}
    for line in open(RUNS / name / ("log.interFoamRhoST" if (RUNS / name / "log.interFoamRhoST").exists() else "log.interFoam"), errors="replace"):
        if line.startswith("deltaT = "):
            cur = {"dt": line.split("=")[1].strip()}
        elif line.startswith("Courant Number mean") and cur is not None:
            cur["co"] = line.strip()
        elif "Phase-1 volume fraction" in line and cur:
            cur["vf"] = line.strip()
        elif line.startswith("ExecutionTime") and cur:
            rows.append(cur)
            cur = {}
            if len(rows) >= n:
                break
    return rows


logs = sorted(p.name for p in (RUNS / "T2_rhoST_check").glob("log.*"))
a_, b_ = steps("T2_rhoST_check", 64), steps("T2_prod", 64)
same = sum(1 for x, y in zip(a_, b_) if x.get("dt") == y.get("dt") and x.get("vf") == y.get("vf"))
out += ["\n### Build check of the rebuilt solver\n",
        "- T2_rhoST_check logs %s; steps read %d and %d; steps with identical deltaT and phase-volume lines: %d" % (logs, len(a_), len(b_), same)]

(HERE / "set6_tests_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
