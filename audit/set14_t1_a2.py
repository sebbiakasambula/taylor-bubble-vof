"""Data audit, set 14: T1med_A2plic, the A2 configuration (isoAdvector with plicRDF and the momentum predictor) at T1 on
the medium mesh, against the T1 baseline V4_medium, the predictor run T1med_mompred and Fluent F_T1, over 5 to 10 s.
Needs audit/m/T1med_A2plic.csv (of_metrics.py) and audit/intact/T1med_A2plic.csv (intact_check.run with tmax 10).
Writes set14_t1_a2_out.md.
"""
import csv
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
LO, HI = 5.0, 10.0


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
    p, cov = np.polyfit(t, x, 1, cov=True) if len(t) > 3 else (np.polyfit(t, x, 1), np.full((2, 2), np.nan))
    return float(p[0]), float(np.sqrt(cov[0, 0])), len(t)


def fr(M):
    w = win(M) & np.isfinite(M["nose"])
    return fit(M["t"][w], M["nose"][w])


def fr_matched(M, B):
    tb = B["t"][win(B) & np.isfinite(B["nose"])]
    idx = [int(np.argmin(np.abs(M["t"] - t))) for t in tb]
    return fit(M["t"][idx], M["nose"][idx]), float(np.max(np.abs(M["t"][idx] - tb)))


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


B, P, A, F = load("V4_medium"), load("T1med_mompred"), load("T1med_A2plic"), load("F_T1")
fb, fp, fa, ff = fr(B), fr(P), fr(A), fr(F)
(fam, _, nam), dt = fr_matched(A, B)
out = ["## Set 14. T1med_A2plic, the A2 configuration at T1 (computed %s)\n" % "14 Nov",
       "Set-up against T1med_mompred: only `application interIsoFoam` and the fvSolution alpha block of T2_A2_plic "
       "(plicRDF). Mesh points and initial alpha identical; endTime 10, writeInterval 0.25 as T1med_mompred.\n",
       "| Fr, 5 to 10 s | Fr | s.e. | n |", "|---|---|---|---|"]
for nm, v in (("V4_medium (baseline)", fb), ("T1med_mompred", fp), ("T1med_A2plic, all rows", fa), ("F_T1 (Fluent)", ff)):
    out.append("| %s | %s | %s | %d |" % (nm, f4(v[0]), f4(v[1]), v[2]))
out.append("| T1med_A2plic, matched to the %d baseline times (max offset %.4f s) | %s | | %d |" % (nam, dt, f4(fam), nam))
out += ["",
        "- A2 against the baseline: %s%% all rows, %s%% matched" % (f4(pct(fa[0], fb[0])), f4(pct(fam, fb[0]))),
        "- A2 against the predictor run: %s%%" % f4(pct(fa[0], fp[0])),
        "- A2 against Fluent: %s%% (baseline %s%%, predictor %s%%)" % (
            f4(pct(fa[0], ff[0])), f4(pct(fb[0], ff[0])), f4(pct(fp[0], ff[0]))),
        "- gap to Fluent as the paper gives it (Fluent / OF - 1): baseline %s%%, predictor %s%%, A2 %s%%" % (
            f4(pct(ff[0], fb[0])), f4(pct(ff[0], fp[0])), f4(pct(ff[0], fa[0]))),
        "- phi (A2 - OF)/(F - OF): %s all rows, %s matched" % (
            f4((fa[0] - fb[0]) / (ff[0] - fb[0])), f4((fam - fb[0]) / (ff[0] - fb[0]))),
        "", "| Quantity, 5 to 10 s | V4_medium | T1med_mompred | T1med_A2plic | A2 vs baseline | F_T1 |",
        "|---|---|---|---|---|---|"]
for k, lab, fn in (("film", "Film (mean)", mean), ("smear", "Smeared gas (mean, %)", mean), ("L_b", "L_b (mean)", mean),
                   ("out_b", "out_b (mean, %)", mean), ("pk_all", "Largest speed (median)", median),
                   ("pk_liq", "Largest liquid speed (median)", median)):
    vb, vp, va = fn(B, k), fn(P, k), fn(A, k)
    try:
        vf = f4(fn(F, k))
    except KeyError:
        vf = "-"
    out.append("| %s | %s | %s | %s | %s%% | %s |" % (lab, f4(vb), f4(vp), f4(va), f4(pct(va, vb)), vf))
g = A["G"]
out += ["", "- gas drift T1med_A2plic, last against t = 0: %s%%" % f4(pct(g[-1], g[0])),
        "- end time %s s, %d snapshots, wall time %s h" % (f4(A["t"][-1]), len(A["t"]), f4(wall("T1med_A2plic")))]
ip = HERE / "intact" / "T1med_A2plic.csv"
if ip.exists():
    R = np.array([[float(x) for x in r] for r in list(csv.reader(open(ip)))[1:]])
    late = R[R[:, 0] >= 1.0]
    out.append("- regions: max %d, max second region %.4f%%, max gas outside the largest after 1 s %.4f%%" % (
        int(R[:, 1].max()), R[:, 3].max(), late[:, 4].max()))
(HERE / "set14_t1_a2_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
