"""Write mirrored structured-grid VTK files of speed and alpha for the fine-mesh split figure (T2, fine mesh):
the fine baseline at 2.6 s and the predictor-only run at 2.6 and 2.8 s, in pipe coordinates (no nose shift).
Same layout as tovtk.py. Output to the pv folder of the work folder."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "audit"))
from of_metrics import RUNS, scal, vec  # noqa: E402

from config import PV as OUT  # noqa: E402
OUT.mkdir(exist_ok=True)


def write(name, x, r, sp, al):
    X = np.unique(np.round(x, 6))
    R = np.unique(np.round(r, 6))
    nx, nr = len(X), len(R)
    assert nx * nr == len(x), (name, nx, nr, len(x))
    ix = np.searchsorted(X, np.round(x, 6))
    ir = np.searchsorted(R, np.round(r, 6))
    S = np.zeros((nr, nx))
    A = np.zeros((nr, nx))
    S[ir, ix] = sp
    A[ir, ix] = al
    Rm = np.r_[-R[::-1], R]
    S = np.vstack([S[::-1], S])
    A = np.vstack([A[::-1], A])
    A[np.abs(A) < 1e-30] = 0.0  # values below single precision (e.g. 1e-47) stop the VTK float reader
    XX, RR = np.meshgrid(X, Rm)
    with open(OUT / (name + ".vtk"), "w") as f:
        f.write("# vtk DataFile Version 3.0\n%s\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d 1\nPOINTS %d float\n"
                % (name, nx, 2 * nr, nx * 2 * nr))
        np.savetxt(f, np.c_[XX.ravel(), RR.ravel(), np.zeros(XX.size)], fmt="%.6g")
        f.write("POINT_DATA %d\nSCALARS speed float 1\nLOOKUP_TABLE default\n" % XX.size)
        np.savetxt(f, np.maximum(S.ravel(), 1e-4), fmt="%.6g")
        f.write("SCALARS alpha float 1\nLOOKUP_TABLE default\n")
        np.savetxt(f, A.ravel(), fmt="%.6g")
    gas = al < 0.5
    print(name, nx, nr, "max gas speed %.2f" % sp[gas].max(), "max liquid speed %.2f" % sp[al > 0.99].max())


for case, t, name in [("T2_C2_fine", "2.6", "split_base_2.6"), ("T2_M2f_mompred_fine", "2.6", "split_pred_2.6"),
                      ("T2_M2f_mompred_fine", "2.8", "split_pred_2.8")]:
    C = vec(RUNS / case / "0/C")
    a = scal(RUNS / case / t / "alpha.water")
    U = vec(RUNS / case / t / "U")
    write(name, C[:, 0], np.round(C[:, 1], 9), np.linalg.norm(U, axis=1), a)
