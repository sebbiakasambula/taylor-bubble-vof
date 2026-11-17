"""R2.1 to R2.8: pre-fixed reading of the T2 density-weighted surface tension test (T2_rhoST), 5 to 8.4 s.
Reads m_T2_rhoST.csv, m_T2_rhoST_all.csv, m_T2_prod.csv, m_F_T2.csv and the t2rho.py arrays; writes t2rho_result.md.
Numbers only. Baseline and Fluent numbers are the fixed values from the plan."""
import csv
from pathlib import Path

import numpy as np

from config import WORK
S = P3 = WORK   # arrays, csv files and the result file are all in the work folder
FR_BASE, FR_FLUENT = 0.218357, 0.263163
FILM_BASE, FILM_FLUENT = 0.134017, 0.142004
FR_MOMPRED = 0.229808
LO, HI = 5.0 - 1e-6, 8.4 + 1e-3


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


def win(M):
    return (M["t"] >= LO) & (M["t"] <= HI)


T = read_csv(P3 / "m_T2_rhoST.csv")
TA = read_csv(P3 / "m_T2_rhoST_all.csv")
B = read_csv(P3 / "m_T2_prod.csv")
F = read_csv(P3 / "m_F_T2.csv")
w, wb, wf = win(T), win(B), win(F)

L = ["# T2_rhoST reading (R2.1 to R2.8), window 5 to 8.4 s", "",
     "Fixed inputs: baseline Fr %.6f, Fluent Fr %.6f, baseline film %.6f, Fluent film %.6f, T2_mompred Fr %.6f. Window rule: %.6f <= t <= %.3f. Standard error = standard error of the least-squares slope." % (
         FR_BASE, FR_FLUENT, FILM_BASE, FILM_FLUENT, FR_MOMPRED, LO, HI), ""]

# self-checks of the reading itself
fb, seb, nb = slope(B["t"][wb], B["nose"][wb])
ff, sef, nf = slope(F["t"][wf], F["nose"][wf])
R = np.load(S / "T2_rhoST_5to84.npy")
mR = (R[:, 0] >= LO) & (R[:, 0] <= HI)
fr_t2rho = float(np.polyfit(R[mR, 0], R[mR, 1], 1)[0])
Rm = np.load(S / "T2_mompred_5to84.npy")
mM = (Rm[:, 0] >= LO) & (Rm[:, 0] <= HI)
fr_mom = float(np.polyfit(Rm[mM, 0], Rm[mM, 1], 1)[0])
phi_mom = (fr_mom - FR_BASE) / (FR_FLUENT - FR_BASE)

fr, se, n = slope(T["t"][w], T["nose"][w])
L += ["## Self-checks", "",
      "- baseline Fr from m_T2_prod.csv: %.6f, n %d, s.e. %.6f (base_5to84.md: 0.218357, 4, 0.001500)" % (fb, nb, seb),
      "- Fluent Fr from m_F_T2.csv: %.6f, n %d, s.e. %.6f (base_5to84.md: 0.263163, 5, 0.001132)" % (ff, nf, sef),
      "- T2_mompred Fr from T2_mompred_5to84.npy: %.6f (expected 0.229808), phi %.4f (expected 0.2556)" % (fr_mom, phi_mom),
      "- t2rho.py on T2_rhoST (T2_rhoST_5to84.npy, n %d): Fr %.9f" % (mR.sum(), fr_t2rho), ""]

# R2.1
d1 = abs(fr - fr_t2rho)
L += ["## R2.1 Fr_test", "",
      "- Fr_test = %.6f (least-squares nose slope, m_T2_rhoST.csv), n = %d, standard error %.6f" % (fr, n, se),
      "- t2rho.py value %.9f; difference %.3e, required <= 1e-9: %s" % (fr_t2rho, d1, "PASS" if d1 <= 1e-9 else "FAIL"),
      "- rows: t %.1f to %.1f" % (T["t"][w].min(), T["t"][w].max()), ""]

# R2.2, R2.3, R2.4
dfr = 100.0 * (fr / FR_BASE - 1.0)
phi = (fr - FR_BASE) / (FR_FLUENT - FR_BASE)
gap = 100.0 * (FR_FLUENT / fr - 1.0)
gap0 = 100.0 * (FR_FLUENT / FR_BASE - 1.0)
L += ["## R2.2 to R2.4", "",
      "- R2.2 dFr = 100 (Fr_test / %.6f - 1) = %.3f %%" % (FR_BASE, dfr),
      "- R2.3 phi = (Fr_test - %.6f) / (%.6f - %.6f) = %.4f (self-test with T2_mompred Fr %.6f: %.4f, expected 0.2556)" % (FR_BASE, FR_FLUENT, FR_BASE, phi, FR_MOMPRED, (FR_MOMPRED - FR_BASE) / (FR_FLUENT - FR_BASE)),
      "- R2.4 remaining gap = 100 (%.6f / Fr_test - 1) = %.2f %%, against %.2f %% before" % (FR_FLUENT, gap, gap0), ""]

# R2.5
f = T["film"][w]
ok = ~np.isnan(f)
mf = float(f[ok].mean())
L += ["## R2.5 film (corrected, wall 0.5 cos 2.5 deg), mean over the window", "",
      "- film = %.6f D, valid rows %d of %d" % (mf, ok.sum(), len(f)),
      "- change against the baseline %.6f: %+.6f D (%+.2f %%)" % (FILM_BASE, mf - FILM_BASE, 100 * (mf / FILM_BASE - 1)),
      "- change against Fluent %.6f: %+.6f D (%+.2f %%)" % (FILM_FLUENT, mf - FILM_FLUENT, 100 * (mf / FILM_FLUENT - 1)),
      "- t2rho.py film (uncorrected, wall 0.5): %.6f; film_raw mean here %.6f" % (float(np.nanmean(R[mR, 2])), float(np.nanmean(T["film_raw"][w]))), ""]

# R2.6
targets = (5.82, 6.66, 7.49, 8.32)
idx = [int(np.argmin(np.abs(T["t"] - tt))) for tt in targets]
fr4, se4, n4 = slope(T["t"][idx], T["nose"][idx])
L += ["## R2.6 refit on the 4 rows nearest 5.82, 6.66, 7.49, 8.32 s", "",
      "- rows used: t = %s" % ", ".join("%.1f" % T["t"][i] for i in idx),
      "- Fr (4 rows) = %.6f, n = %d, standard error %.6f" % (fr4, n4, se4),
      "- compare: Fr (18 rows) = %.6f, standard error %.6f; baseline (4 rows) %.6f, standard error %.6f" % (fr, se, fb, seb),
      "- dFr (4 rows) = %.3f %%, phi (4 rows) = %.4f" % (100 * (fr4 / FR_BASE - 1), (fr4 - FR_BASE) / (FR_FLUENT - FR_BASE)), ""]

# R2.7
G = TA["G"]
drift = 100.0 * (G[-1] / G[0] - 1.0)
L += ["## R2.7 gas drift (m_T2_rhoST_all.csv)", "",
      "- rows %d, t %.1f to %.1f; G first %.9e, G last %.9e" % (len(G), TA["t"][0], TA["t"][-1], G[0], G[-1]),
      "- drift = 100 (G_last / G_first - 1) = %.4f %%; largest |G / G_first - 1| over all rows %.4f %%" % (drift, float(100 * np.max(np.abs(G / G[0] - 1)))), ""]

# R2.8
def means(M, m):
    return {k: float(np.nanmean(M[k][m])) for k in ("L_b", "out_b", "wake_b", "smear")} | {k: float(np.nanmedian(M[k][m])) for k in ("pk_all", "pk_liq", "pk_gas")}


S1, S2, S3 = means(T, w), means(B, wb), means(F, wf)
L += ["## R2.8 field quantities, 5 to 8.4 s rows", "",
      "| quantity | T2_rhoST (n = %d) | T2_prod (n = %d) | Fluent T2 (n = %d) |" % (w.sum(), wb.sum(), wf.sum()), "|---|---|---|---|"]
for k, lab in (("L_b", "mean L_b (D)"), ("out_b", "mean out_b (%)"), ("wake_b", "mean wake_b (%)"), ("smear", "mean smear (%)"),
               ("pk_all", "median pk_all (m/s)"), ("pk_liq", "median pk_liq (m/s)"), ("pk_gas", "median pk_gas (m/s)")):
    L.append("| %s | %.4f | %.4f | %.4f |" % (lab, S1[k], S2[k], S3[k]))
L += ["", "Per-row values, T2_rhoST:", "", "| t | nose | tail_ax | L_b | film | smear % | pk_all |", "|---|---|---|---|---|---|---|"]
for i in np.flatnonzero(w):
    L.append("| %.1f | %.4f | %.4f | %.4f | %.4f | %.3f | %.4f |" % (T["t"][i], T["nose"][i], T["tail_ax"][i], T["L_b"][i], T["film"][i], T["smear"][i], T["pk_all"][i]))
L += ["", "## Decision rule", "",
      "|dFr| = %.3f %%; threshold 2 %%: |dFr| > 2 %% is %s (shift of Fr_test from the baseline is %s)." % (abs(dfr), "TRUE" if abs(dfr) > 2 else "FALSE", "%+.3f %%" % dfr),
      "Combined with the standard errors: Fr_test s.e. %.6f (%.2f %% of Fr_test); baseline s.e. %.6f." % (se, 100 * se / fr, seb), ""]
txt = "\n".join(L) + "\n"
(P3 / "t2rho_result.md").write_text(txt, encoding="utf-8", newline="\n")
print(txt)
