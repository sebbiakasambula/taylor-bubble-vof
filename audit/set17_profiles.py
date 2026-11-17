"""Set 17, part 2. Axial velocity across the film at T2 against Brown's laminar profile (17 Nov).

For each run and snapshot in 5 to 8.4 s, every axial column inside the film window (from 1.5 D below the nose to 0.5 D
above the axis tail, as the film metric) gives the axial velocity of the cells between the wall and the interface. The
theory is eq. (3.5) of Abubakar and Matar with the local film of the column,
    w(r) = (G/4)(R^2 - r^2) - (G a^2/2) ln(R/r),   a = R - delta,   G = drho g / mu_l,
taken downward. Per column: the film flux from the cells against the flux of the theory (both over the liquid part of the
column, r > a), and the speed of the cell nearest the interface against the theory there.

Writes set17_profiles_out.md and profiles.npz (one snapshot per run near 8.3 s, for the figure).
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
sys.path.insert(0, str(HERE))
import tb_post  # noqa: E402
from of_metrics import scal, vec  # noqa: E402
from tb_params import fluid  # noqa: E402

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # shared modules and config.py
from config import FLUENT as FD  # noqa: E402
WALL_OF = 0.5 * math.cos(math.radians(2.5))
f2 = fluid(100, 50, 1000)
G = (f2["rho_l"] - f2["rho_g"]) / f2["mu_l"]
CLEAR = 1.5


def theory(r, R, a):
    return G / 4 * (R * R - r * r) - G * a * a / 2 * np.log(R / r)


def columns(x, y, alpha, ux, wall):
    """Per film-window column: (x, delta, r, w) with r and w of the liquid cells r > a, w downward."""
    nose = tb_post.nose_position(x, y, alpha)
    tail = tb_post.tail_position(x, y, alpha)
    xi, ri = tb_post.interface_radius(x, y, alpha)
    out = []
    for xc, rc in zip(xi, ri):
        if not (tail + 0.5 < xc < nose - CLEAR):
            continue
        m = np.abs(np.round(x, 10) - xc) < 1e-9
        r, w = y[m], -ux[m]
        o = np.argsort(r)
        r, w = r[o], w[o]
        keep = r > rc
        out.append((xc - nose, wall - rc, r[keep], w[keep], r, w, alpha[m][o]))
    return nose, out


def of_snap(name, t):
    case = RUNS / name
    C = vec(case / "0/C")
    d = case / t
    return C[:, 0], np.round(C[:, 1], 9), scal(d / "alpha.water"), vec(d / "U")[:, 0], WALL_OF


def fl_snap(folder, prefix, step):
    dd = np.loadtxt(FD / folder / ("%s-%d" % (prefix, step)), delimiter=",", skiprows=1)
    return dd[:, 1], np.round(dd[:, 2], 9), 1.0 - dd[:, 9], dd[:, 4], 0.5


def stats(cols, R):
    qr, ir = [], []
    for xr, delta, r, w, *_ in cols:
        a = R - delta
        # flux over the liquid cells against the theory over the same cells (midpoint rule on cell centres)
        if len(r) < 3:
            continue
        edges = np.r_[a, 0.5 * (r[1:] + r[:-1]), R]
        dr2 = np.diff(edges ** 2)
        qs = np.sum(w * np.pi * dr2)
        rr = np.linspace(a, R, 4001)
        qt = np.trapezoid(theory(rr, R, a) * 2 * np.pi * rr, rr)
        qr.append(qs / qt)
        ir.append(w[0] / theory(np.array([r[0]]), R, a)[0])
    return np.array(qr), np.array(ir)


def main():
    lines = ["# Set 17, part 2. Film velocity at T2 against Brown's profile", "",
             "G = drho g / mu_l = %.4f. Window: 1.5 D below the nose to 0.5 D above the axis tail." % G, "",
             "| run | snapshots | columns | film flux / theory, median (10th to 90th percentile) | interface-cell speed / theory, median (10th to 90th) | mean film in these columns |",
             "|---|---|---|---|---|---|"]
    keep = {}
    of_runs = [("OF baseline", "T2_prod", ["5.823776", "6.655744", "7.487712", "8.31968"]),
               ("OF predictor", "T2_mompred", ["5.2", "6", "7", "8.2"]),
               ("OF geometric and predictor", "T2_A2_plic", ["5.2", "6", "7", "8.2"])]
    for lab, name, ts in of_runs:
        allq, alli, alld, ncol = [], [], [], 0
        for t in ts:
            x, y, al, ux, wall = of_snap(name, t)
            nose, cols = columns(x, y, al, ux, wall)
            q, i = stats(cols, wall)
            allq += list(q); alli += list(i); alld += [c[1] for c in cols]; ncol += len(cols)
            if t == ts[-1]:
                keep[lab] = (float(t), nose, cols, wall)
        allq, alli = np.array(allq), np.array(alli)
        lines.append("| %s | %s | %d | %.4f (%.4f to %.4f) | %.4f (%.4f to %.4f) | %.4f |" % (
            lab, ", ".join(ts), ncol, np.median(allq), *np.percentile(allq, [10, 90]), np.median(alli),
            *np.percentile(alli, [10, 90]), np.mean(alld)))
    steps = [6160, 7040, 7920, 8800, 9680]
    allq, alli, alld, ncol = [], [], [], 0
    for s in steps:
        x, y, al, ux, wall = fl_snap("T2_case", "T2_Taylor_bubble", s)
        nose, cols = columns(x, y, al, ux, wall)
        q, i = stats(cols, wall)
        allq += list(q); alli += list(i); alld += [c[1] for c in cols]; ncol += len(cols)
        if s == steps[-1]:
            keep["Fluent"] = (8.2, nose, cols, wall)
    allq, alli = np.array(allq), np.array(alli)
    lines.append("| Fluent | steps %s | %d | %.4f (%.4f to %.4f) | %.4f (%.4f to %.4f) | %.4f |" % (
        ", ".join(map(str, steps)), ncol, np.median(allq), *np.percentile(allq, [10, 90]), np.median(alli),
        *np.percentile(alli, [10, 90]), np.mean(alld)))
    # one column in the middle of the window of the last snapshot, for the figure
    lines += ["", "Middle column of the last snapshot of each run (figure):", ""]
    arrays = {}
    for k, (t, nose, cols, wall) in keep.items():
        j = len(cols) // 2
        xr, delta, r, w, rall, wall_w, aall = cols[j]
        a = wall - delta
        lines.append("- %s at t %.4f: column %.4f D behind the nose, film %.4f, %d liquid cells, interface speed %.4f against theory %.4f" % (
            k, t, -xr, delta, len(r), w[0], theory(np.array([r[0]]), wall, a)[0]))
        key = re.sub(r"\W", "_", k)
        arrays[key + "_r"] = rall
        arrays[key + "_w"] = wall_w
        arrays[key + "_alpha"] = aall
        arrays[key + "_meta"] = np.array([t, xr, delta, wall])
    np.savez(HERE / "profiles.npz", **arrays)
    (HERE / "set17_profiles_out.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
