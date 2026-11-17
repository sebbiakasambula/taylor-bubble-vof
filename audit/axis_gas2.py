"""Read-only, part 2 of the fast-gas look (Hint 14). (C) When the fast flow starts: largest Courant number in time
windows over the first 0.3 s. (D) Its spatial structure: share of bubble gas above 2 m/s, and whether the axial
velocity along the axis row and across the bubble changes sign from cell to cell (noise) or varies smoothly (a flow).
Appends to axisgas/axis_gas_out.md."""
import re
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from of_metrics import RUNS, scal, vec  # noqa: E402
from intact_check import merged_index  # noqa: E402
from axis_gas import log_courant  # noqa: E402  (re-runs part A/B printout on import; output file rewritten first)

OUT = HERE / "axisgas" / "axis_gas_out.md"
L = ["\n## C. Onset: largest Courant number in windows of the first 0.3 s\n",
     "| run | " + " | ".join("%.3f-%.3f" % (a, b) for a, b in zip(np.arange(0, 0.3, 0.025), np.arange(0.025, 0.301, 0.025))) + " |",
     "|---|" + "---|" * 12]
for n in ["T2_C2_fine", "T2_M2f_mompred_fine", "T2_A2f_fine"]:
    a = np.array(log_courant(n, 0.31))
    vals = []
    for lo in np.arange(0, 0.3, 0.025):
        s = a[(a[:, 0] > lo) & (a[:, 0] <= lo + 0.025), 1]
        vals.append("%.3f" % s.max() if len(s) else "-")
    L.append("| %s | " % n + " | ".join(vals) + " |")
    # first step above 2 and 5 times the fine baseline's level (0.012)
    for thr in (0.025, 0.05, 0.1):
        hit = a[a[:, 1] > thr]
        L.append("   - %s: first step with Co max > %.3f at t = %s" % (n, thr, ("%.4f s" % hit[0, 0]) if len(hit) else "never"))

L.append("\n## D. Structure of the gas velocity in the main bubble\n")
L.append("Sign changes = share of neighbouring cell pairs whose axial velocity has opposite sign (along the axis row inside "
         "the bubble, and across the bubble in the column of the fastest cell). Jump = mean |difference| between "
         "neighbours over the mean |Ux|. A smooth flow gives few sign changes and a small jump; cell-to-cell noise "
         "gives many.\n")
L.append("| run | t | gas cells in bubble | share > 2 m/s | median gas |U| | 95th pct | axis row: sign changes, jump | "
         "across: sign changes, jump | Ux along the axis row near the peak (10 cells) |")
L.append("|---|---|---|---|---|---|---|---|---|")
for n in ["T2_C2_fine", "T2_M2f_mompred_fine", "T2_A2f_fine"]:
    case = RUNS / n
    C = vec(case / "0/C")
    ux, ix = merged_index(C[:, 0])
    uy, iy = merged_index(np.round(C[:, 1], 9))
    for t in ("0.2", "1", "2"):
        a = scal(case / t / "alpha.water")
        U = vec(case / t / "U")
        A = np.ones((len(ux), len(uy)))
        A[ix, iy] = a
        X = np.zeros_like(A)
        X[ix, iy] = U[:, 0]
        M = np.zeros_like(A)
        M[ix, iy] = np.linalg.norm(U, axis=1)
        lab, nreg = ndimage.label(A < 0.5)
        sizes = ndimage.sum(np.ones_like(A), lab, range(1, nreg + 1))
        body = lab == 1 + int(np.argmax(sizes))
        mg = M[body]
        i0 = np.where(body[:, 0])[0]
        ax = X[i0, 0]

        def stats(v):
            if len(v) < 3:
                return float("nan"), float("nan")
            s = np.mean(np.sign(v[1:]) != np.sign(v[:-1]))
            j = np.mean(np.abs(np.diff(v))) / max(np.mean(np.abs(v)), 1e-12)
            return s, j
        s_ax, j_ax = stats(ax)
        kk = np.unravel_index(np.argmax(np.where(body, M, -1)), M.shape)
        col = X[kk[0], :][body[kk[0], :]]
        s_c, j_c = stats(col)
        p = int(np.argmax(np.abs(ax)))
        win = ax[max(0, p - 5):p + 5]
        L.append("| %s | %s | %d | %.1f%% | %.2f | %.2f | %.2f, %.2f | %.2f, %.2f | %s |" % (
            n, t, body.sum(), 100 * np.mean(mg > 2), np.median(mg), np.percentile(mg, 95), s_ax, j_ax, s_c, j_c,
            " ".join("%+.1f" % v for v in win)))

with open(OUT, "a", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print("\n".join(L))
