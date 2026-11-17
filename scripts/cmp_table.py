"""B8: comparison tables from the per-snapshot metric CSVs (OpenFOAM m_*_prod / m_V4_medium / m_VA_square, Fluent m_F_*).

Inputs and outputs (all in the work folder, $P2_WORK): cmp_table.csv / .md, t1mesh_table.csv / .md, wake_t3.csv / .md, base_5to84.csv / .md.
Window rows: 5 - 1e-6 <= t <= 20.001 (base table: <= 8.4 + 1e-3). Fr = least-squares slope of nose against t.
Percentages are computed from unrounded values, 100 * (a / b - 1). Nothing in metrics.py or m_*.csv is changed.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import csv
import ctypes
import sys
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import percentileofscore

sys.stdout.reconfigure(encoding="utf-8")
try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)  # BELOW_NORMAL
except Exception:
    pass
from config import WORK

FILES = {"OF_T1": "m_V4_medium.csv", "OF_T1sq": "m_VA_square.csv", "OF_T2": "m_T2_prod.csv", "OF_T3": "m_T3_prod.csv",
         "OF_T3all": "m_T3_prod_all.csv", "F_T1": "m_F_T1.csv", "F_T1sq": "m_F_T1sq.csv", "F_T2": "m_F_T2.csv", "F_T3": "m_F_T3.csv"}
VIANA = {"T1": 0.338404, "T2": 0.250410, "T3": 0.277025}
LLEW = {"T1": 0.061508, "T2": 0.139409}
LO, HI = 5 - 1e-6, 20.001


def load(f):
    with open(WORK / f, newline="") as fh:
        rd = list(csv.DictReader(fh))
    return {k: (np.array([r[k] for r in rd]) if k == "pk_region" else np.array([float(r[k]) for r in rd])) for k in rd[0]}


D = {k: load(f) for k, f in FILES.items()}


def window(M, hi=HI, lo=LO):
    return (M["t"] >= lo) & (M["t"] <= hi)


def pct(a, b):
    return 100.0 * (a / b - 1.0)


def fr(M, w):
    """Least-squares slope of nose against t over the rows w, its standard error, and the number of rows."""
    ok = w & np.isfinite(M["nose"]); t, y = M["t"][ok], M["nose"][ok]; n = int(ok.sum())
    p = np.polyfit(t, y, 1)
    se = math_sqrt(np.sum((y - np.polyval(p, t)) ** 2) / (n - 2) / np.sum((t - t.mean()) ** 2)) if n > 2 else float("nan")
    return float(p[0]), se, n


def math_sqrt(x):
    return float(np.sqrt(x))


def fin(M, key, w):
    v = M[key][w]
    return v[np.isfinite(v)]


def mean_n(M, key, w):
    v = fin(M, key, w)
    return float(v.mean()), len(v)


def med_n(M, key, w):
    v = fin(M, key, w)
    return float(np.median(v)), len(v)


# ---------- formatting ----------
def f(v, d):
    return "n/a" if v is None or not np.isfinite(v) else "%.*f" % (d, v)


def sg(v, d=1):
    if v is None or not np.isfinite(v): return "n/a"
    s = "%+.*f" % (d, v)
    return s if float(s) != 0 else "%.*f" % (d, 0)


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def write_csv(name, header, recs):
    with open(WORK / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(header)
        for r in recs: w.writerow([repr(float(x)) if isinstance(x, (float, np.floating)) else x for x in r])


# =====================================================================================================================
# 1. cmp_table: three cases, OpenFOAM and Fluent, medium mesh
# =====================================================================================================================
recs = []   # id, case, code, statistic, value, n
R = {}      # (id, statistic, case, code) -> (value, n)


def add(i, case, code, stat, v, n):
    recs.append((i, case, code, stat, float(v), int(n))); R[(i, stat, case, code)] = (float(v), int(n))


for case in ("T1", "T2", "T3"):
    for code in ("OF", "F"):
        M = D["%s_%s" % (code, case)]; w = window(M)
        s, se, n = fr(M, w)
        add("C1", case, code, "Fr", s, n)
        add("C2", case, code, "Fr_vs_Viana_pct", pct(s, VIANA[case]), n)
        v, n = mean_n(M, "film", w); add("C4", case, code, "film_mean", v, n)
        if case in LLEW: add("C5", case, code, "film_vs_Llewellin_pct", pct(v, LLEW[case]), n)
        for stat, key, i in (("L_b_mean", "L_b", "C7"), ("L_b_min", "L_b", "C7"), ("L_b_max", "L_b", "C7")):
            x = fin(M, key, w); add(i, case, code, stat, {"mean": x.mean, "min": x.min, "max": x.max}[stat.split("_")[-1]](), len(x))
        for stat, key, i in (("out_b_mean", "out_b", "C9"), ("wake_b_mean", "wake_b", "C10"), ("smear_mean", "smear", "C11")):
            v, n = mean_n(M, key, w); add(i, case, code, stat, v, n)
        for stat, key, i in (("pk_all_median", "pk_all", "C12"), ("pk_liq_median", "pk_liq", "C13"), ("pk_gas_median", "pk_gas", "C14")):
            v, n = med_n(M, key, w); add(i, case, code, stat, v, n)
        add("C15", case, code, "n_window_rows", int(w.sum()), int(w.sum()))
    add("C3", case, "F_vs_OF", "Fr_F_vs_OF_pct", pct(R[("C1", "Fr", case, "F")][0], R[("C1", "Fr", case, "OF")][0]), R[("C1", "Fr", case, "F")][1])
    add("C6", case, "F_vs_OF", "film_F_vs_OF_pct", pct(R[("C4", "film_mean", case, "F")][0], R[("C4", "film_mean", case, "OF")][0]), R[("C4", "film_mean", case, "F")][1])
write_csv("cmp_table.csv", ["id", "case", "code", "statistic", "value", "n"], recs)


def cells(i, stat, fmt, kind="both", with_n=False):
    out = []
    for case in ("T1", "T2", "T3"):
        for code in ("OF", "F"):
            if kind == "ratio":
                if code == "OF": out.append("–"); continue
                v, n = R[(i, stat, case, "F_vs_OF")]
            else:
                key = (i, stat, case, code)
                if key not in R: out.append("n/a"); continue
                v, n = R[key]
            out.append(fmt(v) + (" (%d)" % n if with_n else ""))
    return out


f4 = lambda v: f(v, 4); f2 = lambda v: f(v, 2); f1 = lambda v: f(v, 1); f0 = lambda v: f(v, 0)
head = ["Quantity", "T1 OF", "T1 F", "T2 OF", "T2 F", "T3 OF", "T3 F"]
rows = [
    ["C1 Fr (nose speed, m/s per D/s units of the run)"] + cells("C1", "Fr", f4),
    ["C2 Fr vs Viana (%)"] + cells("C2", "Fr_vs_Viana_pct", sg),
    ["C3 Fr, F vs OF (%)"] + cells("C3", "Fr_F_vs_OF_pct", sg, "ratio"),
    ["C4 film (D), valid rows in brackets"] + cells("C4", "film_mean", f4, with_n=True),
    ["C5 film vs Llewellin (%)"] + cells("C5", "film_vs_Llewellin_pct", sg),
    ["C6 film, F vs OF (%)"] + cells("C6", "film_F_vs_OF_pct", sg, "ratio"),
    ["C7 L_b mean (D)"] + cells("C7", "L_b_mean", f2), ["C7 L_b min (D)"] + cells("C7", "L_b_min", f2), ["C7 L_b max (D)"] + cells("C7", "L_b_max", f2),
    ["C9 gas outside the body, mean (% of gas)"] + cells("C9", "out_b_mean", f1),
    ["C10 gas in the 2 D wake, mean (% of gas)"] + cells("C10", "wake_b_mean", f1),
    ["C11 smeared gas, mean (% of gas)"] + cells("C11", "smear_mean", f1),
    ["C12 median peak speed, all cells (m/s)"] + cells("C12", "pk_all_median", f2),
    ["C13 median peak speed, liquid cells (m/s)"] + cells("C13", "pk_liq_median", f2),
    ["C14 median peak speed, gas cells (m/s)"] + cells("C14", "pk_gas_median", f2),
    ["C15 window rows"] + cells("C15", "n_window_rows", f0),
]
rows[0][0] = "C1 Fr"
NOTE = ("Window 5 to 20 s (5 - 1e-6 <= t <= 20.001). Medium mesh: OF = V4_medium, T2_prod, T3_prod; F = Fluent T1, T2, T3. "
        "Fr = least-squares slope of the nose against t. Viana Fr 0.338404, 0.250410, 0.277025; Llewellin film 0.061508, 0.139409 (none for T3). "
        "Percentages from unrounded values, F vs OF = 100 (F / OF - 1).")
(WORK / "cmp_table.md").write_text("# cmp_table\n\n" + table(head, rows) + "\n\n" + NOTE + "\n", encoding="utf-8")

# =====================================================================================================================
# 2. t1mesh_table: the four T1 runs
# =====================================================================================================================
T1 = [("OF stretched", "OF_T1"), ("OF square", "OF_T1sq"), ("F stretched", "F_T1"), ("F square", "F_T1sq")]
recs2 = []; Q = {}


def add2(i, run, stat, v, n):
    recs2.append((i, run, stat, float(v), int(n))); Q[(i, stat, run)] = (float(v), int(n))


for run, key in T1:
    M = D[key]; w = window(M)
    s, se, n = fr(M, w); add2("M1", run, "Fr", s, n)
    v, n = mean_n(M, "film", w); add2("M2", run, "film_mean", v, n)
    x = fin(M, "L_b", w); add2("M3", run, "L_b_mean", x.mean(), len(x)); add2("M3", run, "L_b_min", x.min(), len(x)); add2("M3", run, "L_b_max", x.max(), len(x))
    for k in ("out_b", "wake_b"):
        x = fin(M, k, w); add2("M4", run, k + "_mean", x.mean(), len(x)); add2("M4", run, k + "_max", x.max(), len(x))
    for k in ("pk_all", "pk_liq", "pk_gas"):
        v, n = med_n(M, k, w); add2("M5", run, k + "_median", v, n)
    v, n = med_n(M, "u_core", w); add2("M6", run, "u_core_median", v, n)
    v, n = med_n(M, "fast_core", w); add2("M7", run, "fast_core_median", v, n)
    cnt = Counter(M["pk_region"][w])
    for k in ("bubble", "wake", "ahead", "beside", "none"): add2("M8", run, "pk_region_" + k, cnt.get(k, 0), int(w.sum()))
    v, n = mean_n(M, "gap05", w); add2("M9", run, "gap05_mean", v, n)
    x = fin(M, "skirt_b", w); add2("M10", run, "skirt_b_mean", x.mean(), len(x)); add2("M10", run, "skirt_b_gt_0.1_share_pct", 100.0 * (x > 0.1).sum() / len(x), len(x))
    add2("M10", run, "skirt_b_gt_0.1_count", int((x > 0.1).sum()), len(x))
# square vs stretched per code
for code, a, b in (("OF", "OF square", "OF stretched"), ("F", "F square", "F stretched")):
    add2("M1", a, "Fr_sq_vs_str_pct", pct(Q[("M1", "Fr", a)][0], Q[("M1", "Fr", b)][0]), Q[("M1", "Fr", a)][1])
    add2("M3", a, "L_b_mean_sq_minus_str_D", Q[("M3", "L_b_mean", a)][0] - Q[("M3", "L_b_mean", b)][0], Q[("M3", "L_b_mean", a)][1])
    add2("M3", a, "L_b_mean_sq_vs_str_pct", pct(Q[("M3", "L_b_mean", a)][0], Q[("M3", "L_b_mean", b)][0]), Q[("M3", "L_b_mean", a)][1])
    for k in ("pk_all", "pk_liq", "pk_gas"):
        add2("M5", a, k + "_sq_vs_str_pct", pct(Q[("M5", k + "_median", a)][0], Q[("M5", k + "_median", b)][0]), Q[("M5", k + "_median", a)][1])

# M11 figure snapshots
SNAPS = [("OF stretched", "OF_T1", "t", 19.700352), ("OF square", "OF_T1sq", "t", 19.700352), ("F stretched", "F_T1", "step", 21645), ("F square", "F_T1sq", "step", 8190)]
M11 = ["L_b", "out_b", "wake_b", "pk_all", "pk_liq", "pk_gas", "u_core", "fast_core", "skirt_b"]
snapinfo = {}
for run, key, how, val in SNAPS:
    M = D[key]; w = window(M)
    idx = np.flatnonzero(np.abs(M["t"] - val) < 1e-6) if how == "t" else np.flatnonzero(M["step"] == val)
    assert len(idx) == 1 and w[idx[0]], (run, idx)
    i = idx[0]; snapinfo[run] = (M["t"][i], int(M["step"][i]) if "step" in M else None)
    add2("M11", run, "t", M["t"][i], int(w.sum()))
    for k in M11:
        ref = fin(M, k, w); v = M[k][i]
        add2("M11", run, k, v, len(ref)); add2("M11", run, k + "_pctrank", percentileofscore(ref, v, kind="rank"), len(ref))
write_csv("t1mesh_table.csv", ["id", "run", "statistic", "value", "n"], recs2)

runs = [r for r, _ in T1]
q = lambda i, s, r: Q[(i, s, r)][0]
qn = lambda i, s, r: Q[(i, s, r)][1]


def rowq(label, i, stat, fmt, only=None, with_n=False):
    out = [label]
    for r in runs:
        if only and r not in only: out.append("–"); continue
        out.append(fmt(q(i, stat, r)) + (" (%d)" % qn(i, stat, r) if with_n else ""))
    return out


SQ = ("OF square", "F square")
rows2 = [
    rowq("M1 Fr", "M1", "Fr", f4), rowq("M1 Fr, square vs stretched (%)", "M1", "Fr_sq_vs_str_pct", sg, SQ),
    rowq("M2 film (D), valid rows in brackets", "M2", "film_mean", f4, with_n=True),
    rowq("M3 L_b mean (D)", "M3", "L_b_mean", f2), rowq("M3 L_b min (D)", "M3", "L_b_min", f2), rowq("M3 L_b max (D)", "M3", "L_b_max", f2),
    rowq("M3 L_b mean, square minus stretched (D)", "M3", "L_b_mean_sq_minus_str_D", lambda v: sg(v, 2), SQ),
    rowq("M3 L_b mean, square vs stretched (%)", "M3", "L_b_mean_sq_vs_str_pct", sg, SQ),
    rowq("M4 out_b mean (% of gas)", "M4", "out_b_mean", f1), rowq("M4 out_b max (% of gas)", "M4", "out_b_max", f1),
    rowq("M4 wake_b mean (% of gas)", "M4", "wake_b_mean", f1), rowq("M4 wake_b max (% of gas)", "M4", "wake_b_max", f1),
    rowq("M5 median pk_all (m/s)", "M5", "pk_all_median", f2), rowq("M5 median pk_liq (m/s)", "M5", "pk_liq_median", f2), rowq("M5 median pk_gas (m/s)", "M5", "pk_gas_median", f2),
    rowq("M5 pk_all, square vs stretched (%)", "M5", "pk_all_sq_vs_str_pct", sg, SQ), rowq("M5 pk_liq, square vs stretched (%)", "M5", "pk_liq_sq_vs_str_pct", sg, SQ),
    rowq("M5 pk_gas, square vs stretched (%)", "M5", "pk_gas_sq_vs_str_pct", sg, SQ),
    rowq("M6 median u_core (m/s)", "M6", "u_core_median", f2, with_n=True), rowq("M7 median fast_core (% of core gas)", "M7", "fast_core_median", f1, with_n=True),
    ["M8 pk_region counts"] + [", ".join("%s %d" % (k, q("M8", "pk_region_" + k, r)) for k in ("bubble", "wake", "ahead", "beside", "none") if q("M8", "pk_region_" + k, r) > 0) for r in runs],
    rowq("M9 mean gap05 (D)", "M9", "gap05_mean", f4),
    rowq("M10 mean skirt_b (D)", "M10", "skirt_b_mean", f2),
    ["M10 share of rows with skirt_b > 0.1 (%, rows)"] + ["%s (%d of %d)" % (f1(q("M10", "skirt_b_gt_0.1_share_pct", r)), q("M10", "skirt_b_gt_0.1_count", r), qn("M10", "skirt_b_gt_0.1_count", r)) for r in runs],
]
fmt11 = {"L_b": f2, "out_b": f1, "wake_b": f1, "pk_all": f2, "pk_liq": f2, "pk_gas": f2, "u_core": f2, "fast_core": f1, "skirt_b": f2}
unit11 = {"L_b": "(D)", "out_b": "(% of gas)", "wake_b": "(% of gas)", "pk_all": "(m/s)", "pk_liq": "(m/s)", "pk_gas": "(m/s)", "u_core": "(m/s)", "fast_core": "(%)", "skirt_b": "(D)"}
rows3 = [["time (s)"] + [f(q("M11", "t", r), 4) for r in runs],
         ["code snapshot"] + ["t = 19.700352 s", "t = 19.700352 s", "step 21645", "step 8190"]]
for k in M11:
    rows3.append(["%s %s" % (k, unit11[k])] + ["%s (%d)" % (fmt11[k](q("M11", k, r)), round(q("M11", k + "_pctrank", r))) for r in runs])
NOTE2 = ("Window 5 to 20 s. OF stretched = V4_medium, OF square = VA_square, F stretched = Fluent T1_Case, F square = Fluent T1 square (with completion). "
         "Fr = least-squares slope of the nose against t. Square vs stretched = 100 (square / stretched - 1) from unrounded values. "
         "u_core and fast_core are medians over the rows where they are defined (count in brackets). "
         "M11: value of the snapshot used in the figures, with its percentile rank (scipy percentileofscore, kind 'rank', 0 to 100) among the window rows of the same run, in brackets.")
(WORK / "t1mesh_table.md").write_text("# t1mesh_table\n\n## M1 to M10\n\n" + table(["Quantity"] + runs, rows2) + "\n\n## M11 figure snapshots, value (percentile rank)\n\n" + table(["Quantity"] + runs, rows3) + "\n\n" + NOTE2 + "\n", encoding="utf-8")

# =====================================================================================================================
# 3. wake_t3
# =====================================================================================================================
recs3 = []; Wk = {}


def add3(i, item, code, stat, v, n):
    recs3.append((i, item, code, stat, float(v), int(n))); Wk[(i, item, stat, code)] = (float(v), int(n))


MO, MF = D["OF_T3"], D["F_T3"]
iO = np.flatnonzero(np.abs(MO["t"] - 20.305026) < 1e-6); iF = np.flatnonzero(MF["step"] == 78000)
assert len(iO) == 1 and len(iF) == 1
W1 = ["wake_ax", "wake_b", "out_ax", "out_b", "wk_nfast_ax", "wk_fastgas_ax", "wk_nfast_b", "wk_fastgas_b"]
for code, M, i in (("OF", MO, iO[0]), ("F", MF, iF[0])):
    add3("W1", "snapshot", code, "t", M["t"][i], 1)
    for k in W1: add3("W1", "snapshot", code, k, M[k][i], 1)
for code, M in (("OF", MO), ("F", MF)):
    w = window(M); v, n = mean_n(M, "wake_b", w); add3("W2", "window_5to20", code, "wake_b_mean", v, n)
    ok = w & (M["wk_nfast_b"] >= 10); x = M["wk_fastgas_b"][ok]; x = x[np.isfinite(x)]
    add3("W3", "window_5to20", code, "wk_fastgas_b_median_nfast_ge_10", np.median(x) if len(x) else float("nan"), len(x))
A = D["OF_T3all"]
for tt in (5, 20, 40):
    i = int(np.argmin(np.abs(A["t"] - tt)))
    add3("W4", "nearest_t_%d" % tt, "OF", "t", A["t"][i], 1); add3("W4", "nearest_t_%d" % tt, "OF", "L_b", A["L_b"][i], 1); add3("W4", "nearest_t_%d" % tt, "OF", "smear", A["smear"][i], 1)
for lab, hi in (("5to20", 20.001), ("5to40", 40.001)):
    w = window(A, hi); x = fin(A, "L_b", w)
    add3("W4", "window_" + lab, "OF", "L_b_min", x.min(), len(x)); add3("W4", "window_" + lab, "OF", "L_b_max", x.max(), len(x))
w = window(A, 40.001); x = fin(A, "smear", w)
add3("W4", "window_5to40", "OF", "smear_min", x.min(), len(x)); add3("W4", "window_5to40", "OF", "smear_max", x.max(), len(x))
write_csv("wake_t3.csv", ["id", "item", "code", "statistic", "value", "n"], recs3)

g = lambda i, item, s, c: Wk[(i, item, s, c)][0]
w1rows = []
for k in W1:
    d = 0 if k.startswith("wk_n") else 2
    w1rows.append([k + (" (% of gas)" if k.startswith(("wake", "out")) or "fastgas" in k else " (cells)")] + [f(g("W1", "snapshot", k, c), d) for c in ("OF", "F")])
nW = lambda i, item, s, c: Wk[(i, item, s, c)][1]
w2 = [["mean wake_b, 5 to 20 s (% of gas)"] + ["%s (%d rows)" % (f(g("W2", "window_5to20", "wake_b_mean", c), 3), nW("W2", "window_5to20", "wake_b_mean", c)) for c in ("OF", "F")]]
w3 = [["median wk_fastgas_b over rows with wk_nfast_b >= 10 (%)"] + ["%s (%d rows)" % (f(g("W3", "window_5to20", "wk_fastgas_b_median_nfast_ge_10", c), 1), nW("W3", "window_5to20", "wk_fastgas_b_median_nfast_ge_10", c)) for c in ("OF", "F")]]
w4a = [[("t = %g s" % tt)] + ["%.3f" % g("W4", "nearest_t_%d" % tt, "t", "OF"), f(g("W4", "nearest_t_%d" % tt, "L_b", "OF"), 2), f(g("W4", "nearest_t_%d" % tt, "smear", "OF"), 1)] for tt in (5, 20, 40)]
w4b = [["L_b min (D)", f(g("W4", "window_5to20", "L_b_min", "OF"), 2), f(g("W4", "window_5to40", "L_b_min", "OF"), 2)],
       ["L_b max (D)", f(g("W4", "window_5to20", "L_b_max", "OF"), 2), f(g("W4", "window_5to40", "L_b_max", "OF"), 2)],
       ["smear min (%)", "–", f(g("W4", "window_5to40", "smear_min", "OF"), 1)], ["smear max (%)", "–", f(g("W4", "window_5to40", "smear_max", "OF"), 1)],
       ["rows", str(nW("W4", "window_5to20", "L_b_min", "OF")), str(nW("W4", "window_5to40", "L_b_min", "OF"))]]
NOTE3 = ("OF = m_T3_prod.csv (W1 to W3) and m_T3_prod_all.csv (W4); F = Fluent T3. W1: OF at t = 20.305026 s, F at step 78000 (t = %.4f s). "
         "Windows: 5 - 1e-6 <= t <= 20.001 (W2, W3, and 5 to 20 in W4), <= 40.001 (5 to 40 in W4); the OF row at 20.305 s is outside the 5 to 20 s window." % g("W1", "snapshot", "t", "F"))
(WORK / "wake_t3.md").write_text("# wake_t3\n\n## W1 snapshot\n\n" + table(["Quantity", "OF 20.305 s", "F step 78000"], w1rows) +
                                 "\n\n## W2 mean wake_b over the window\n\n" + table(["Quantity", "OF", "F"], w2) +
                                 "\n\n## W3 fast cells in the wake containing gas\n\n" + table(["Quantity", "OF", "F"], w3) +
                                 "\n\n## W4 OpenFOAM T3 over the long run (m_T3_prod_all.csv)\n\n" + table(["Row nearest", "t (s)", "L_b (D)", "smear (%)"], w4a) +
                                 "\n\n" + table(["Quantity", "5 to 20 s", "5 to 40 s"], w4b) + "\n\n" + NOTE3 + "\n", encoding="utf-8")

# =====================================================================================================================
# 4. base_5to84
# =====================================================================================================================
recs4 = []; Bs = {}
HI84 = 8.4 + 1e-3
for case in ("T2", "T3"):
    for code in ("OF", "F"):
        M = D["%s_%s" % (code, case)]; w = window(M, HI84)
        s, se, n = fr(M, w)
        vals = [("Fr", s, n), ("Fr_se", se, n)]
        for stat, key in (("film_mean", "film"), ("out_b_mean", "out_b"), ("wake_b_mean", "wake_b"), ("smear_mean", "smear")):
            v, nn = mean_n(M, key, w); vals.append((stat, v, nn))
        for stat, v, nn in vals:
            recs4.append(("B", case, code, stat, float(v), int(nn))); Bs[(case, code, stat)] = (float(v), int(nn))
write_csv("base_5to84.csv", ["id", "case", "code", "statistic", "value", "n"], recs4)
b = lambda c, k, s: Bs[(c, k, s)]
rows4 = []
for case in ("T2", "T3"):
    for code in ("OF", "F"):
        rows4.append(["%s %s" % (case, code), f(b(case, code, "Fr")[0], 6), str(b(case, code, "Fr")[1]), f(b(case, code, "Fr_se")[0], 6), f(b(case, code, "film_mean")[0], 6),
                      f(b(case, code, "out_b_mean")[0], 4), f(b(case, code, "wake_b_mean")[0], 4), f(b(case, code, "smear_mean")[0], 2)])
NOTE4 = "Window 5 - 1e-6 <= t <= 8.4 + 1e-3. Fr = least-squares slope of the nose against t with n rows and the slope standard error. Film, out_b, wake_b (% of gas) and smear (% of gas) are means over the window rows."
(WORK / "base_5to84.md").write_text("# base_5to84\n\n" + table(["Run", "Fr", "n", "Fr s.e.", "film (D)", "out_b (%)", "wake_b (%)", "smear (%)"], rows4) + "\n\n" + NOTE4 + "\n", encoding="utf-8")

for n_ in ("cmp_table.md", "t1mesh_table.md", "wake_t3.md", "base_5to84.md"):
    print((WORK / n_).read_text(encoding="utf-8"))
