"""Read-only look at the fast gas next to the axis (Hint 14). Part A: Courant number history from the logs, to see when
the fast flow starts. Part B: the fastest gas cells in the reconstructed fields, with their row from the axis, position
behind the nose and velocity components. Writes axisgas/axis_gas_out.md."""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from of_metrics import RUNS, scal, vec  # noqa: E402
from intact_check import merged_index  # noqa: E402

OUT = HERE / "axisgas"
OUT.mkdir(exist_ok=True)
RUNS_A = ["T2_C2_fine", "T2_M2f_mompred_fine", "T2_A2f_fine", "T2_mompred", "T2_prod", "T2_C1_coarse",
          "T2_M1c_mompred_coarse"]
lines = []


def log_courant(name, tmax=1.0):
    log = RUNS / name / "log.interFoam"
    if not log.exists():
        log = RUNS / name / ("log.interIsoFoam")
    co, out = None, []
    with open(log, errors="replace") as f:
        for ln in f:
            if ln.startswith("Courant Number mean"):
                m = re.search(r"max:\s*([0-9.eE+-]+)", ln)
                co = float(m[1]) if m else None
            elif ln.startswith("Time = ") and co is not None:
                t = float(ln.split("=")[1])
                out.append((t, co))
                co = None
                if t > tmax:
                    break
    return out


def main():
    lines.append("## A. Largest cell Courant number from the logs (computed at the start of each step)\n")
    lines.append("| run | first step | max to 0.01 s | max to 0.05 s | max to 0.2 s | max to 1 s | first t with Co > 1 |")
    lines.append("|---|---|---|---|---|---|---|")
    for n in RUNS_A:
        try:
            h = log_courant(n)
        except FileNotFoundError:
            lines.append("| %s | no log |" % n)
            continue
        a = np.array(h)

        def mx(t):
            s = a[a[:, 0] <= t + 1e-12, 1]
            return s.max() if len(s) else float("nan")
        big = a[a[:, 1] > 1.0]
        lines.append("| %s | %.4f (t %.2e) | %.3f | %.3f | %.3f | %.3f | %s |" % (
            n, a[0, 1], a[0, 0], mx(0.01), mx(0.05), mx(0.2), mx(1.0), ("%.4f s" % big[0, 0]) if len(big) else "none"))

    lines.append("\n## B. Fastest gas cells (alpha < 0.5) in the reconstructed fields\n")
    lines.append("Row 0 is the cell row next to the axis. x - nose in D (negative = behind the nose). Ux axial (up is +), "
                 "Ur radial.\n")
    for n in ["T2_C2_fine", "T2_M2f_mompred_fine", "T2_A2f_fine", "T2_mompred"]:
        case = RUNS / n
        C = vec(case / "0/C")
        ux, ix = merged_index(C[:, 0])
        uy, iy = merged_index(np.round(C[:, 1], 9))
        r = np.hypot(C[:, 1], C[:, 2])
        times = [t for t in ("0.2", "1", "2", "2.6") if (case / t / "U").exists()]
        lines.append("### %s (%d rows from axis to wall, first row centre r %.4f D, its height %.4f D)" % (
            n, len(uy), uy[0], 2 * uy[0]))
        lines.append("| t | nose x | max gas |U| | row | x - nose | Ux | Ur | cells with gas |U| > 2 m/s, by row (0,1,2,3,4+) | "
                     "max liquid |U| |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for t in times:
            a = scal(case / t / "alpha.water")
            U = vec(case / t / "U")
            mag = np.linalg.norm(U, axis=1)
            gas = a < 0.5
            # nose: highest x of a gas cell in row 0
            g0 = gas & (iy == 0)
            nose = ux[ix[g0]].max() if g0.any() else float("nan")
            k = np.where(gas)[0][np.argmax(mag[gas])]
            urad = (U[k, 1] * C[k, 1] + U[k, 2] * C[k, 2]) / max(r[k], 1e-12)
            fast = gas & (mag > 2.0)
            rows = [int((fast & (iy == j)).sum()) for j in range(4)] + [int((fast & (iy >= 4)).sum())]
            lines.append("| %s | %.3f | %.2f | %d | %.3f | %+.2f | %+.2f | %s | %.2f |" % (
                t, nose, mag[k], iy[k], ux[ix[k]] - nose, U[k, 0], urad, rows, mag[~gas].max()))
        lines.append("")

    text = "\n".join(lines) + "\n"
    (OUT / "axis_gas_out.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
