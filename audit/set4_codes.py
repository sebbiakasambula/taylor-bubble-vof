"""Data audit, set 4: OpenFOAM against Fluent on the medium mesh at T1 to T3 (V4_medium, T2_prod, T3_prod; Fluent T1,
T2, T3), the runs table, and the Fluent T1 run continued with 40 iterations per step.

Uses the regenerated cmp/cmp_table.csv and cmp/wake_t3.csv (cmp_table.py on the recomputed per-snapshot tables), the
recomputed tables in audit/m, the solver logs, and the residual histories inside the Fluent .dat.h5 files.
Writes set4_codes_out.md.
"""
import csv
import glob
import math
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
# run_new.py sits in this folder
import run_new  # noqa: E402

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import FLUENT as FD  # noqa: E402
LO = 5.0 - 1e-6
VIANA = {"T1": 0.338404, "T2": 0.250410, "T3": 0.277025}
LLEW = {"T1": 0.061508, "T2": 0.139409}
GCI_MED, GCI_FINE = 2.3043, 1.0577  # set 2


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


def fit(t, y):
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    n = len(t)
    return float(p[0]), math.sqrt(np.sum(r ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)), n


def frw(M, lo, hi):
    t = M["t"]
    w = (t >= lo) & (t <= hi + 1e-3) & np.isfinite(M["nose"])
    return fit(t[w], M["nose"][w])


out = ["## Set 4. OpenFOAM against Fluent on the medium mesh, recomputed 14 Nov\n",
       "cmp_table.csv and wake_t3.csv regenerated from the recomputed per-snapshot tables (identical to the stored ones, see set 3). "
       "Window 5 to 20 s unless stated.\n", "### cmp_table (regenerated), 4 decimal places\n",
       "| id | case | code | statistic | value | n |", "|---|---|---|---|---|---|"]
C = list(csv.DictReader(open(HERE / "cmp" / "cmp_table.csv")))
for r in C:
    try:
        v = f4(float(r["value"]))
    except ValueError:
        v = r["value"]
    out.append("| %s | %s | %s | %s | %s | %s |" % (r["id"], r["case"], r["code"], r["statistic"], v, r["n"]))
out += ["\n### wake_t3 (regenerated), 4 decimal places\n", "| id | item | code | statistic | value | n |", "|---|---|---|---|---|---|"]
for r in csv.DictReader(open(HERE / "cmp" / "wake_t3.csv")):
    try:
        v = f4(float(r["value"]))
    except ValueError:
        v = r["value"]
    out.append("| %s | %s | %s | %s | %s | %s |" % tuple(r.get(k, "") for k in ("id", "item", "code", "statistic")) + (v, r.get("n", "")) if False else
               "| %s | %s | %s | %s | %s | %s |" % (r.get("id", ""), r.get("item", ""), r.get("code", ""), r.get("statistic", ""), v, r.get("n", "")))

val = {(r["case"], r["code"], r["statistic"]): float(r["value"]) for r in C if r["value"] not in ("", "nan")}
out += ["\n### Derived here\n"]
for case in ("T1", "T2", "T3"):
    o, f = val[(case, "OF", "Fr")], val[(case, "F", "Fr")]
    d = 100 * (f / o - 1)
    out.append("- %s: Fr OF %s, F %s, F against OF %s%% (%s times the medium GCI %s%%); against Viana %s: OF %s%% (|dFr| %s), F %s%% (|dFr| %s)" % (
        case, f4(o), f4(f), f4(d), f4(abs(d) / GCI_MED), f4(GCI_MED), f4(VIANA[case]), f4(100 * (o / VIANA[case] - 1)),
        f4(abs(o - VIANA[case])), f4(100 * (f / VIANA[case] - 1)), f4(abs(f - VIANA[case]))))
    fo, ff = val[(case, "OF", "film_mean")], val[(case, "F", "film_mean")]
    s = "  film OF %s, F %s, F against OF %s%%" % (f4(fo), f4(ff), f4(100 * (ff / fo - 1)))
    if case in LLEW:
        s += "; against Llewellin %s: OF %s%%, F %s%%" % (f4(LLEW[case]), f4(100 * (fo / LLEW[case] - 1)), f4(100 * (ff / LLEW[case] - 1)))
    out.append(s)
    lo_, lf = val[(case, "OF", "L_b_mean")], val[(case, "F", "L_b_mean")]
    out.append("  bubble length OF %s, F %s, F against OF %s%%" % (f4(lo_), f4(lf), f4(100 * (lf / lo_ - 1))))
out.append("- Mesh shape at T1 against the medium GCI: 3.5784 / %s = %s; codes at T1: %s" % (
    f4(GCI_MED), f4(3.5784 / GCI_MED), f4(abs(100 * (val[("T1", "F", "Fr")] / val[("T1", "OF", "Fr")] - 1)) / GCI_MED)))

M1f = load("V4_fine")
w = (M1f["t"] >= LO) & (M1f["t"] <= 20.001) & np.isfinite(M1f["film"])
frf = frw(M1f, LO, 1e9)[0]
out.append("- T1 fine mesh: Fr 5 s to end %s, against Viana %s%%; film 5 to 20 s %s, against Llewellin %s%%" % (
    f4(frf), f4(100 * (frf / VIANA["T1"] - 1)), f4(M1f["film"][w].mean()), f4(100 * (M1f["film"][w].mean() / LLEW["T1"] - 1))))

# long runs
T2 = load("T2_prod")
a = frw(T2, LO, 1e9)
out.append("- T2_prod nose 5 s to end (t last %s): Fr %s (n %d)" % (f4(T2["t"][-1]), f4(a[0]), a[2]))
T3 = load("T3_prod")
b = frw(T3, LO, 40.0)
k40 = np.argmin(np.abs(T3["t"] - 40.0))
out.append("- T3_prod nose 5 to 40 s: Fr %s (n %d), against Viana %s%%; nose at t = %s: %s D" % (
    f4(b[0]), b[2], f4(100 * (b[0] / VIANA["T3"] - 1)), f4(T3["t"][k40]), f4(T3["nose"][k40])))
for tt in (5.3, 39.9):
    k = np.argmin(np.abs(T3["t"] - tt))
    out.append("- T3_prod at t = %s: L_b %s D, smear %s%%" % (f4(T3["t"][k]), f4(T3["L_b"][k]), f4(T3["smear"][k])))
for tt in (5.0, 40.0):
    k = np.argmin(np.abs(T3["t"] - tt))
    out.append("- T3_prod at t = %s: smear %s%%" % (f4(T3["t"][k]), f4(T3["smear"][k])))
w5 = (T3["t"] >= LO) & (T3["t"] <= 20.001)
out.append("- T3_prod L_b over 5 to 20 s: %s to %s D (n %d)" % (f4(np.nanmin(T3["L_b"][w5])), f4(np.nanmax(T3["L_b"][w5])), int(w5.sum())))


# centroid of T2_prod from the raw fields
def scal(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(", s)
    return np.fromstring(s[m.end():s.index(")", m.end())], sep=" ", count=int(m[1]))


def vec(f):
    s = Path(f).read_text()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*(\d+)\s*\(", s)
    n = int(m[1])
    body = s[m.end():]
    body = body[:body.index("\n)")].replace("(", " ").replace(")", " ")
    return np.fromstring(body, sep=" ", count=3 * n).reshape(n, 3)


case = RUNS / "T2_prod"
Cx = vec(case / "0/C")[:, 0]
V = scal(case / "0/V")
tc, xc = [], []
for d in sorted([d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.]+", d.name)], key=lambda d: float(d.name)):
    t = float(d.name)
    if t < LO:
        continue
    g = (1 - scal(d / "alpha.water")) * V
    tc.append(t)
    xc.append(float(np.sum(g * Cx) / np.sum(g)))
cf = fit(np.array(tc), np.array(xc))
out.append("- T2_prod gas centroid 5 s to end: Fr %s (n %d)" % (f4(cf[0]), cf[2]))


# runs table
def wall_time(name):
    ct = re.findall(r"ClockTime = ([0-9.]+) s", (RUNS / name / "log.interFoam").read_text(errors="replace"))
    return float(ct[-1]) / 3600


def maxdt(name):
    return float(re.search(r"maxDeltaT\s+([0-9.eE+-]+)", (RUNS / name / "system/controlDict").read_text())[1])


out += ["\n### Runs table\n", "| case | OF end (s) | OF maxDeltaT | OF wall time (h) | OF largest gas drift (%) | Fluent last snapshot (s) | Fluent final data (s) | Fluent step cap | cap against OF (%) | Fluent largest gas drift (%) |",
        "|---|---|---|---|---|---|---|---|---|---|"]
for case, of, ff, fold in (("T1", "V4_medium", "F_T1", "T1_Case"), ("T2", "T2_prod", "F_T2", "T2_case"), ("T3", "T3_prod", "F_T3", "T3_Case")):
    Mo, Mf = load(of), load(ff)
    go = 100 * np.max(np.abs(Mo["G"] / Mo["G"][0] - 1))
    gf = 100 * np.max(np.abs(Mf["G"] / Mf["G"][0] - 1))
    anc = run_new.dat_anchors(FD / fold)
    outr = run_new.out_rows(FD / fold)
    last = max(list(anc.values()) + [max(outr.values())])
    s = sorted(outr)[-200:]
    cap = float(np.median(np.diff([outr[k] for k in s]) / np.diff(s)))
    out.append("| %s | %s | %.6g | %s | %s | %s | %s | %.6g | %s | %.3g |" % (case, f4(Mo["t"][-1]), maxdt(of), f4(wall_time(of)), f4(go), f4(Mf["t"][-1]),
                                                                       f4(last), cap, f4(100 * (cap / maxdt(of) - 1)), gf))


# Fluent continuity residual at the end of each step (numcheck1/g_resid.py method)
def endvals(path, nit):
    f = h5py.File(path, "r")
    g = f["results/residuals/phase-1/continuity"]
    d = g["data"][:]
    it = g["iterations"][:]
    k = len(it) - 1
    while k > 0 and it[k - 1] == it[k] - 1:
        k -= 1
    d, it = d[k:], it[k:]
    ends = [i for i in range(len(it)) if int(it[i]) % nit == 0]
    return d[ends, 0] / d[ends, 1][0], int(it[-1])


out += ["\n### Fluent continuity residual at the end of each step (scaled, from the residual history in the .dat.h5)\n"]
for lab, p, nit in (("T1, 20 iterations", "T1_Case/T1_Taylor_bubble.dat.h5", 20),
                    ("T1, 40 iterations (step 24000)", "40_iterations/T1_Taylor_bubble-10-24000.dat.h5", 40),
                    ("T1 square", "T1_square_w_completion/Taylor_bubble_square.dat.h5", 20),
                    ("T2", "T2_case/T2_Taylor_bubble-9-20000.dat.h5", 20),
                    ("T3", "T3_Case/T3_Taylor_bubble-9-85000.dat.h5", 20)):
    try:
        s, itl = endvals(FD / p, nit)
        out.append("- %s: %d step ends, median %.4g, mean %.4g, min %.3g, max %.3g (last iteration counter %d)" % (
            lab, len(s), np.median(s), np.mean(s), s.min(), s.max(), itl))
    except Exception as e:  # noqa: BLE001
        out.append("- %s: not readable (%s)" % (lab, e))

# 40 iterations
F40, F1 = load("F_T1_40it"), load("F_T1")
a40 = fit(F40["t"], F40["nose"])
a20 = frw(F1, LO, 20.0)
wf = (F1["t"] >= LO) & (F1["t"] <= 20.001) & np.isfinite(F1["film"])
out += ["\n### Fluent T1 with 40 iterations per step\n",
        "- %d snapshots, t %s to %s; Fr %s (s.e. %s, %s%% of Fr); 20 iterations over 5 to 20 s: %s; difference %s%%" % (
            len(F40["t"]), f4(F40["t"][0]), f4(F40["t"][-1]), f4(a40[0]), "%.6f" % a40[1], f4(100 * a40[1] / a40[0]), f4(a20[0]),
            f4(100 * (a40[0] / a20[0] - 1))),
        "- film: 40 iterations %s (n %d), 20 iterations 5 to 20 s %s; change %s%%" % (
            f4(np.nanmean(F40["film"])), int(np.isfinite(F40["film"]).sum()), f4(F1["film"][wf].mean()),
            f4(100 * (np.nanmean(F40["film"]) / F1["film"][wf].mean() - 1)))]

(HERE / "set4_codes_out.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out[out.index("### Derived here\n") if "### Derived here\n" in out else 0:]))
