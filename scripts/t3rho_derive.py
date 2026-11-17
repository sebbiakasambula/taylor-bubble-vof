"""R3.1 to R3.9: pre-fixed reading of a T3 density-weighted surface tension test, window 5 to 8.4 s.

python t3rho_derive.py [CASE] [--all-t1 49] [--no-run]

CASE defaults to T3_rhoST. Steps: run_of.py (subprocess) over 5 to 8.4 s with --clear 1.5 into m_<CASE>_5to84.csv and over
all times (0.1 to --all-t1) into m_<CASE>_all.csv; read the T3 baselines from base_5to84.csv, m_T3_prod.csv and m_F_T3.csv;
write t3rho_result.md (CASE T3_rhoST) or t3rho_result_<CASE>.md. With --all-t1 other than 49 the all-times file is named
m_<CASE>_all_to<t1>.csv so that an existing full file is not overwritten. Numbers only.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

from config import WORK
P3 = WORK   # csv files and the result file are all in the work folder
T0, T1, CLEAR = 5.0, 8.4, 1.5
LO, HI = T0 - 1e-6, T1 + 1e-3


def read_csv(p):
    with open(p, newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        rows = list(rd)
    M = {}
    for i, k in enumerate(head):
        col = [r[i] for r in rows]
        try:
            M[k] = np.array([float(v) for v in col])
        except ValueError:
            M[k] = np.array(col)
    return M


def slope(t, y):
    """Least-squares slope, its standard error and n."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    p = np.polyfit(t, y, 1)
    res = y - np.polyval(p, t)
    se = float(np.sqrt(np.sum(res ** 2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))) if len(t) > 2 else float("nan")
    return float(p[0]), se, len(t)


def baselines():
    """T3 rows of base_5to84.csv as {(code, statistic): value}."""
    out = {}
    with open(P3 / "base_5to84.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["case"] == "T3":
                out[(r["code"], r["statistic"])] = float(r["value"])
    return out


def win(M):
    return (M["t"] >= LO) & (M["t"] <= HI)


def paths(case, all_t1):
    f_win = P3 / ("m_%s_5to84.csv" % case)
    f_all = P3 / (("m_%s_all.csv" % case) if all_t1 == 49 else ("m_%s_all_to%g.csv" % (case, all_t1)))
    return f_win, f_all


def run_of(case, all_t1):
    f_win, f_all = paths(case, all_t1)
    for t0, t1, out in ((T0, T1, f_win), (0.1, all_t1, f_all)):
        cmd = [sys.executable, str(Path(__file__).resolve().with_name("run_of.py")), "--case", case, "--t0", str(t0), "--t1", str(t1), "--clear", str(CLEAR), "--out", str(out)]
        print(" ".join(cmd), flush=True)
        r = subprocess.run(cmd, cwd=str(P3))
        if r.returncode != 0:
            sys.exit("run_of.py failed for %s" % out.name)


def compute(case, all_t1=49):
    """All R3 numbers as a dict, from the CSV files (which must exist)."""
    f_win, f_all = paths(case, all_t1)
    base = baselines()
    T, TA = read_csv(f_win), read_csv(f_all)
    B, F = read_csv(P3 / "m_T3_prod.csv"), read_csv(P3 / "m_F_T3.csv")
    w, wb, wf = win(T), win(B), win(F)
    R = dict(case=case, files=(f_win.name, f_all.name), base=base, base_n=base_n(), w=w, wb=wb, wf=wf, T=T, TA=TA, B=B, F=F)

    # R3.1 to R3.3
    fr, se, n = slope(T["t"][w], T["nose"][w])
    fr0, fr_f = base[("OF", "Fr")], base[("F", "Fr")]
    R.update(fr=fr, se=se, n=n, dfr=100 * (fr / fr0 - 1) + 0.0, phi=(fr - fr0) / (fr_f - fr0) + 0.0, gap=100 * (fr_f / fr - 1), gap0=100 * (fr_f / fr0 - 1))

    # R3.4
    f = T["film"][w]
    ok = ~np.isnan(f)
    R.update(film=float(f[ok].mean()), film_n=int(ok.sum()), film_rows=int(w.sum()))

    # R3.5 and R3.8
    def stats(M, m):
        d = {k: float(np.nanmean(M[k][m])) for k in ("out_b", "wake_b", "smear", "L_b")}
        d.update(L_b_min=float(np.nanmin(M["L_b"][m])), L_b_max=float(np.nanmax(M["L_b"][m])))
        d.update({k + "_med": float(np.nanmedian(M[k][m])) for k in ("pk_all", "pk_liq", "pk_gas")})
        return d
    R.update(s_test=stats(T, w), s_prod=stats(B, wb), s_fl=stats(F, wf))

    # R3.6
    G = TA["G"]
    R.update(drift=100 * (G[-1] / G[0] - 1), drift_max=float(100 * np.max(np.abs(G / G[0] - 1))), G0=float(G[0]), G1=float(G[-1]),
             nall=len(G), t_all=(float(TA["t"][0]), float(TA["t"][-1])))

    # R3.7: rows nearest to the T3_prod snapshot times inside the window
    tt = B["t"][wb]
    idx = [int(np.argmin(np.abs(T["t"] - v))) for v in tt]
    fr7, se7, n7 = slope(T["t"][idx], T["nose"][idx])
    R.update(fr7=fr7, se7=se7, n7=n7, idx7=idx, tt=tt, dfr7=100 * (fr7 / fr0 - 1) + 0.0, phi7=(fr7 - fr0) / (fr_f - fr0) + 0.0)

    # R3.9: last row of the window
    def last(M, m):
        i = np.flatnonzero(m)[-1]
        return float(M["t"][i]), float(M["wake_b"][i]), float(M["wk_fastgas_b"][i]), float(M["wk_nfast_b"][i])
    R.update(last_test=last(T, w), last_prod=last(B, wb), last_fl=last(F, wf))
    return R


def write_md(R, path):
    b, fr0, fr_f = R["base"], R["base"][("OF", "Fr")], R["base"][("F", "Fr")]
    T, w = R["T"], R["w"]
    st, sp, sf = R["s_test"], R["s_prod"], R["s_fl"]
    L = ["# %s reading (R3.1 to R3.9), window %g to %g s" % (R["case"], T0, T1), "",
         "Files: %s, %s. Baselines from base_5to84.csv: T3 OpenFOAM Fr %.6f, Fluent Fr %.6f, films %.6f and %.6f. Window rule: %.6f <= t <= %.3f. Standard error = standard error of the least-squares slope." % (
             R["files"][0], R["files"][1], fr0, fr_f, b[("OF", "film_mean")], b[("F", "film_mean")], LO, HI), "",
         "## R3.1 Fr_test", "",
         "- Fr_test = %.6f (least-squares nose slope), n = %d, standard error %.6f (%.2f %% of Fr_test); rows t %.4f to %.4f" % (
             R["fr"], R["n"], R["se"], 100 * R["se"] / R["fr"], T["t"][w].min(), T["t"][w].max()),
         "- baseline: Fr %.6f, n %d, standard error %.6f; Fluent: Fr %.6f, n %d, standard error %.6f" % (
             fr0, int(R["base_n"]["OF"]), b[("OF", "Fr_se")], fr_f, int(R["base_n"]["F"]), b[("F", "Fr_se")]), "",
         "## R3.2 to R3.3", "",
         "- R3.2 dFr = 100 (Fr_test / %.6f - 1) = %.3f %%; |dFr| > 2 %%: %s" % (fr0, R["dfr"], "TRUE" if abs(R["dfr"]) > 2 else "FALSE"),
         "- R3.3 phi = (Fr_test - %.6f) / (%.6f - %.6f) = %.4f" % (fr0, fr_f, fr0, R["phi"]),
         "- remaining gap = 100 (%.6f / Fr_test - 1) = %.2f %%, against %.2f %% before" % (fr_f, R["gap"], R["gap0"]), "",
         "## R3.4 film (corrected, wall 0.5 cos 2.5 deg), mean over the window", "",
         "- film = %.6f D, valid rows %d of %d" % (R["film"], R["film_n"], R["film_rows"]),
         "- change against OpenFOAM %.6f: %+.6f D (%+.2f %%)" % (b[("OF", "film_mean")], R["film"] - b[("OF", "film_mean")], 100 * (R["film"] / b[("OF", "film_mean")] - 1)),
         "- change against Fluent %.6f: %+.6f D (%+.2f %%)" % (b[("F", "film_mean")], R["film"] - b[("F", "film_mean")], 100 * (R["film"] / b[("F", "film_mean")] - 1)), "",
         "## R3.5 and R3.8 field quantities, window rows", "",
         "| quantity | %s test (n = %d) | T3_prod baseline (n = %d) | Fluent T3 (n = %d) |" % (R["case"], w.sum(), R["wb"].sum(), R["wf"].sum()), "|---|---|---|---|"]
    for k, lab in (("out_b", "mean out_b (%)"), ("wake_b", "mean wake_b (%)"), ("smear", "mean smear (%)"), ("L_b", "mean L_b (D)"),
                   ("L_b_min", "min L_b (D)"), ("L_b_max", "max L_b (D)"), ("pk_all_med", "median pk_all (m/s)"),
                   ("pk_liq_med", "median pk_liq (m/s)"), ("pk_gas_med", "median pk_gas (m/s)")):
        L.append("| %s | %.4f | %.4f | %.4f |" % (lab, st[k], sp[k], sf[k]))
    L += ["", "base_5to84.csv, OpenFOAM T3: out_b %.4f, wake_b %.4f, smear %.4f; Fluent T3: out_b %.4f, wake_b %.4f, smear %.4f" % (
        b[("OF", "out_b_mean")], b[("OF", "wake_b_mean")], b[("OF", "smear_mean")], b[("F", "out_b_mean")], b[("F", "wake_b_mean")], b[("F", "smear_mean")]), "",
          "## R3.6 gas drift (%s)" % R["files"][1], "",
          "- rows %d, t %.4f to %.4f; G first %.9e, G last %.9e" % (R["nall"], R["t_all"][0], R["t_all"][1], R["G0"], R["G1"]),
          "- drift = 100 (G_last / G_first - 1) = %.4f %%; largest |G / G_first - 1| over all rows %.4f %%" % (R["drift"], R["drift_max"]), "",
          "## R3.7 refit on the rows nearest to the T3_prod snapshot times in the window", "",
          "- T3_prod times: %s" % ", ".join("%.4f" % v for v in R["tt"]),
          "- rows used: t = %s (%d distinct)" % (", ".join("%.4f" % T["t"][i] for i in R["idx7"]), len(set(R["idx7"]))),
          "- Fr = %.6f, n = %d, standard error %.6f; dFr = %.3f %%, phi = %.4f" % (R["fr7"], R["n7"], R["se7"], R["dfr7"], R["phi7"]), "",
          "## R3.9 wake at the last window row", "",
          "| run | t | wake_b (%) | wk_fastgas_b (%) | fast cells (U > 0.3) |", "|---|---|---|---|---|"]
    for lab, v in ((R["case"] + " test", R["last_test"]), ("T3_prod baseline", R["last_prod"]), ("Fluent T3", R["last_fl"])):
        L.append("| %s | %.4f | %.4f | %.3f | %d |" % (lab, v[0], v[1], v[2], v[3]))
    L += ["", "Per-row values:", "", "| t | nose | tail_b | L_b | film | smear % | out_b % | wake_b % | pk_all |", "|---|---|---|---|---|---|---|---|---|"]
    for i in np.flatnonzero(w):
        L.append("| %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f | %.3f | %.4f |" % (T["t"][i], T["nose"][i], T["tail_b"][i], T["L_b"][i], T["film"][i], T["smear"][i], T["out_b"][i], T["wake_b"][i], T["pk_all"][i]))
    Path(path).write_text("\n".join(L) + "\n", encoding="utf-8", newline="\n")
    return "\n".join(L) + "\n"


def base_n():
    n = {}
    with open(P3 / "base_5to84.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["case"] == "T3" and r["statistic"] == "Fr":
                n[r["code"]] = int(r["n"])
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="T3_rhoST")
    ap.add_argument("--all-t1", type=float, default=49)
    ap.add_argument("--no-run", action="store_true")
    a = ap.parse_args()
    if not a.no_run:
        run_of(a.case, a.all_t1)
    R = compute(a.case, a.all_t1)
    out = P3 / ("t3rho_result.md" if a.case == "T3_rhoST" else "t3rho_result_%s.md" % a.case)
    print(write_md(R, out))
    print("written", out)


if __name__ == "__main__":
    main()
