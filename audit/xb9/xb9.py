"""B9: independent check of the T3_rhoST test run against the T3_prod baseline (window 5 to 8.4 s).

Own code: plain text parsing of the OpenFOAM ASCII fields, scipy.ndimage for the connected regions.
No import of metrics.py, run_of.py, t3rho_derive.py, tb_post.py or snap.py.
Adapted from xb5.py (nose rule, film rule, slope and its standard error are the same code).
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import ctypes
import math
import re
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.stdout.reconfigure(encoding="utf-8")
try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)  # BELOW_NORMAL
except Exception:
    pass
HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))  # shared modules and config.py
from config import RUNS  # noqa: E402
T_LO, T_HI = 5.0 - 1e-6, 8.4 + 1e-3
CLEAR = 1.5
MIN_COLS = 3
WALL_SHIFT = 0.5 - 0.5 * math.cos(math.radians(2.5))
EXPECT = {"Fr": 0.320527, "Fr_se": 0.010495, "n": 5, "film": 0.077342, "out_b": 3.9306, "wake_b": 2.3187, "smear": 4.34}
S4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])      # 4-connectivity


def f6(v):
    return "nan" if v is None or not np.isfinite(v) else "%.6g" % v


def read_field(path, kind):
    """internalField of an ASCII OpenFOAM field: 'scalar' gives (N,), 'vector' gives (N, 3)."""
    txt = Path(path).read_text()
    i = txt.index("internalField")
    m = re.compile(r"internalField\s+nonuniform\s+List<%s>\s*(\d+)\s*\(" % kind).search(txt, i)
    if m is None:
        u = re.compile(r"internalField\s+uniform\s+([^;]+);").search(txt, i)
        raise ValueError("%s: no nonuniform list (%s)" % (path, u.group(0) if u else "unknown format"))
    n = int(m.group(1))
    body = txt[m.end():]
    body = body[:body.index("\n)")]
    if kind == "vector":
        v = np.array(body.replace("(", " ").replace(")", " ").split(), float).reshape(-1, 3)
    else:
        v = np.array(body.split(), float)
    if len(v) != n:
        raise ValueError("%s: expected %d values, read %d" % (path, n, len(v)))
    return v


def time_folders(case, lo=T_LO, hi=T_HI):
    out = []
    for d in Path(case).iterdir():
        if d.is_dir():
            try:
                t = float(d.name)
            except ValueError:
                continue
            if lo <= t <= hi:
                out.append((t, d))
    return sorted(out, key=lambda p: p[0])


def axis_row(x, y):
    """Indices of the axis cells (smallest y), sorted by x."""
    idx = np.flatnonzero(y <= y.min() + 1e-9)
    return idx[np.argsort(x[idx])]


def nose_and_tail(x, y, a):
    """Same rule as xb5.py: nose = interpolated alpha = 0.5 crossing after the highest-x axis cell with alpha < 0.5;
    tail = interpolated crossing before the last consecutive run of axis cells with alpha < 0.5 (None if the run reaches x = 0)."""
    ia = axis_row(x, y)
    xa, aa = x[ia], a[ia]
    gas = np.flatnonzero(aa < 0.5)
    if len(gas) == 0:
        return None, None, 0
    k = int(gas[-1])                                  # highest-x axis cell with alpha < 0.5
    if k + 1 < len(xa):
        nose = xa[k] + (0.5 - aa[k]) * (xa[k + 1] - xa[k]) / (aa[k + 1] - aa[k])
    else:
        nose = xa[k]
    j = k
    while j > 0 and aa[j - 1] < 0.5:                  # walk down to the last consecutive cell with alpha < 0.5
        j -= 1
    tail = None
    if j > 0:
        tail = xa[j - 1] + (0.5 - aa[j - 1]) * (xa[j] - xa[j - 1]) / (aa[j] - aa[j - 1])
    return float(nose), (None if tail is None else float(tail)), len(ia)


def slope_se(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float); n = len(t)
    tm, ym = t.mean(), y.mean()
    sxx = np.sum((t - tm) ** 2)
    b = np.sum((t - tm) * (y - ym)) / sxx
    a0 = ym - b * tm
    rss = np.sum((y - a0 - b * t) ** 2)
    return float(b), float(math.sqrt(rss / (n - 2) / sxx)) if n > 2 else float("nan"), n


def film_of(xr, y, a, nose, tail):
    """Mean of 0.5 - r(interface) over the columns with tail + 0.5 < x < nose - CLEAR.
    Returns (mean, number of columns used); mean is nan when fewer than MIN_COLS columns are usable."""
    order = np.lexsort((y, xr))                       # by column, then by radius
    xs, ys, as_ = xr[order], y[order], a[order]
    cuts = np.flatnonzero(np.diff(xs) != 0) + 1
    vals = []
    for cx, cy, ca in zip(np.split(xs, cuts), np.split(ys, cuts), np.split(as_, cuts)):
        x0 = cx[0]
        if not (tail + 0.5 < x0 < nose - CLEAR):
            continue
        cr = np.flatnonzero((ca[:-1] < 0.5) & (ca[1:] >= 0.5))
        if len(cr) == 0:
            continue
        j = int(cr[-1])                               # last crossing, going outward
        r = cy[j] + (0.5 - ca[j]) * (cy[j + 1] - cy[j]) / (ca[j + 1] - ca[j])
        vals.append(0.5 - r)
    if len(vals) < MIN_COLS:
        return float("nan"), len(vals)
    return float(np.mean(vals)), len(vals)


class Grid:
    """Structured x-r layout of the flat cell list: idx2d[i, j] is the flat index of axial column i, radial row j."""

    def __init__(self, case):
        C = read_field(RUNS / case / "0" / "C", "vector")
        self.V = read_field(RUNS / case / "0" / "V", "scalar")
        self.x, self.y = C[:, 0], C[:, 1]
        self.xr = np.round(self.x, 6)
        self.ux, ix = np.unique(self.xr, return_inverse=True)
        self.uy, iy = np.unique(np.round(self.y, 9), return_inverse=True)
        self.nx, self.nr = len(self.ux), len(self.uy)
        if self.nx * self.nr != len(self.x):
            raise ValueError("%s: %d cells but %d x %d grid" % (case, len(self.x), self.nx, self.nr))
        self.idx2d = np.full((self.nx, self.nr), -1, dtype=int)
        self.idx2d[ix, iy] = np.arange(len(self.x))
        if (self.idx2d < 0).any() or len(np.unique(self.idx2d)) != len(self.x):
            raise ValueError("%s: cells do not form a structured grid" % case)
        self.V2 = self.V[self.idx2d]
        self.X2 = np.broadcast_to(self.ux[:, None], (self.nx, self.nr))
        self.ncell = len(self.x)


def body_metrics(g, a, nose):
    """Body region = 4-connected region of cells with alpha < 0.5 that contains the nose cell
    (highest-x axis-row cell with alpha < 0.5). All other cells are 'detached'."""
    A2 = a[g.idx2d]
    gas2 = A2 < 0.5
    kax = np.flatnonzero(gas2[:, 0])
    lab, nlab = ndimage.label(gas2, structure=S4)
    if len(kax) == 0:                                 # no axis cell with alpha < 0.5: no nose, no body region
        gasvol = (1.0 - A2) * g.V2
        G = float(gasvol.sum())
        nan = float("nan")
        return dict(tail_b=nan, L_b=nan, out_b=nan, wake_b=nan, smear=100.0 * float(gasvol[A2 > 0.5].sum()) / G, G=G,
                    ncomp=int(nlab), nose_is_biggest=False, body_gas_pct=nan, body_cells=0, ax_gap=0, nose_cell_x=nan, ax_body_last_x=nan, det_gas_pct_all=nan)
    k = int(kax[-1])                                  # nose cell, same as the axis rule used for the nose
    body = lab == lab[k, 0]
    ax_body = np.flatnonzero(body[:, 0])
    tail_b = float(g.ux[ax_body[0]])
    gasvol = (1.0 - A2) * g.V2
    G = float(gasvol.sum())
    det = ~body
    out_b = 100.0 * float(gasvol[det & (g.X2 < tail_b - 0.05)].sum()) / G
    wake_b = 100.0 * float(gasvol[det & (g.X2 > tail_b - 2.0) & (g.X2 < tail_b - 0.05)].sum()) / G
    smear = 100.0 * float(gasvol[A2 > 0.5].sum()) / G
    # diagnostics
    comp_gas = ndimage.sum(gasvol * gas2, lab, index=np.arange(1, nlab + 1))
    biggest = int(np.argmax(comp_gas)) + 1
    ax_gap = int(ax_body[-1] - ax_body[0] + 1 - len(ax_body))   # axis cells between tail_b and nose cell missing from the body
    return dict(tail_b=tail_b, L_b=nose - tail_b, out_b=out_b, wake_b=wake_b, smear=smear, G=G,
                ncomp=int(nlab), nose_is_biggest=bool(lab[k, 0] == biggest), body_gas_pct=100.0 * float(gasvol[body].sum()) / G,
                body_cells=int(body.sum()), ax_gap=ax_gap, nose_cell_x=float(g.ux[k]), ax_body_last_x=float(g.ux[ax_body[-1]]),
                det_gas_pct_all=100.0 * float(gasvol[det].sum()) / G)


def run(case, lo, hi):
    g = Grid(case)
    rows = []
    for t, d in time_folders(RUNS / case, lo, hi):
        a = read_field(d / "alpha.water", "scalar")
        nose, tail, naxis = nose_and_tail(g.x, g.y, a)
        f, nc = film_of(g.xr, g.y, a, nose, tail) if (nose is not None and tail is not None) else (float("nan"), 0)
        r = dict(t=t, name=d.name, nose=(float("nan") if nose is None else nose), tail_xb5=tail, film=f, ncols=nc)
        r.update(body_metrics(g, a, nose))
        rows.append(r)
    return g, rows


def summarise(rows):
    ok = [r for r in rows if np.isfinite(r["nose"])]            # folders with at least one axis cell with alpha < 0.5
    t = np.array([r["t"] for r in ok]); nose = np.array([r["nose"] for r in ok])
    fr, se, n = slope_se(t, nose)
    films = np.array([r["film"] for r in rows]); good = np.isfinite(films)
    film_raw = float(films[good].mean()) if good.any() else float("nan")
    Lb = np.array([r["L_b"] for r in ok])
    return dict(Fr=fr, se=se, n=n, nfold=len(rows), film_raw=film_raw, film=film_raw - WALL_SHIFT, nfilm=int(good.sum()),
                smear=float(np.mean([r["smear"] for r in rows])), out_b=float(np.mean([r["out_b"] for r in ok])),
                wake_b=float(np.mean([r["wake_b"] for r in ok])), Lb_mean=float(Lb.mean()), Lb_min=float(Lb.min()), Lb_max=float(Lb.max()),
                Lb_min_t=float(t[Lb.argmin()]), Lb_max_t=float(t[Lb.argmax()]),
                no_nose=[r["name"] for r in rows if not np.isfinite(r["nose"])])


def table(rows):
    L = ["| t (s) | nose | tail_b | L_b | film (0.5 - r) | film columns | smear (%) | out_b (%) | wake_b (%) |", "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| %s | %s | %s | %s | %s | %d | %s | %s | %s |" % (r["name"], f6(r["nose"]), f6(r["tail_b"]), f6(r["L_b"]), f6(r["film"]), r["ncols"],
                                                              f6(r["smear"]), f6(r["out_b"]), f6(r["wake_b"])))
    return L


def diag_table(rows):
    L = ["| t (s) | tail_x (xb5 rule) | body components | nose comp. is largest | body gas share (%) | gap cells on axis in body | nose cell x | last axis body cell x |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| %s | %s | %d | %s | %s | %d | %s | %s |" % (r["name"], f6(r["tail_xb5"]), r["ncomp"], "yes" if r["nose_is_biggest"] else "NO",
                                                          f6(r["body_gas_pct"]), r["ax_gap"], f6(r["nose_cell_x"]), f6(r["ax_body_last_x"])))
    return L


gs, rs = run("T3_rhoST", T_LO, T_HI)
gp, rp = run("T3_prod", T_LO, T_HI)
ss, sp = summarise(rs), summarise(rp)
shift = 100.0 * (ss["Fr"] / sp["Fr"] - 1.0)

# total gas at every time folder 0.2 ... 8.4 of T3_rhoST (and folder 0 for reference)
gas_rows = []
for t, d in time_folders(RUNS / "T3_rhoST", -1.0, 8.4 + 1e-3):
    a = read_field(d / "alpha.water", "scalar")
    gas_rows.append((t, d.name, float(((1.0 - a)[gs.idx2d] * gs.V2).sum())))
g_win = [r for r in gas_rows if 0.2 - 1e-6 <= r[0]]
G_first, G_last = g_win[0][2], g_win[-1][2]
drift = 100.0 * (G_last / G_first - 1.0)
Gs = np.array([r[2] for r in g_win]); dG = 100.0 * np.diff(Gs) / Gs[:-1]

L = ["# B9 check: T3_rhoST (test) and T3_prod (baseline), independent code (`xb9.py`)", "",
     "Own parser and own formulas; no use of metrics.py, run_of.py, t3rho_derive.py, tb_post.py or snap.py. Window 5 - 1e-6 <= t <= 8.4 + 1e-3. Numbers to 6 significant figures.", "",
     "## Grid", "",
     "- T3_rhoST: %d cells, %d axial columns x %d radial rows (x from %s to %s, axis row y = %s, outermost row y = %s)." % (gs.ncell, gs.nx, gs.nr, f6(gs.ux[0]), f6(gs.ux[-1]), f6(gs.uy[0]), f6(gs.uy[-1])),
     "- T3_prod: %d cells, %d axial columns x %d radial rows (same x and y lists: %s)." % (gp.ncell, gp.nx, gp.nr, "yes" if (gp.nx == gs.nx and gp.nr == gs.nr and np.allclose(gp.ux, gs.ux, atol=1e-9) and np.allclose(gp.uy, gs.uy, atol=1e-9)) else "NO"),
     "- Wall in the wedge mid-plane: 0.5 cos(2.5 deg) = %s; film wall shift 0.5 - 0.5 cos(2.5 deg) = %s." % (f6(0.5 * math.cos(math.radians(2.5))), f6(WALL_SHIFT)), "",
     "## Rules used", "",
     "- nose: interpolated alpha = 0.5 crossing after the highest-x axis cell with alpha < 0.5 (same rule as xb5.py). Fr = least-squares slope of nose against t.",
     "- film: columns with tail + 0.5 < x < nose - %g (tail from the xb5.py axis walk), last outward crossing of alpha = 0.5 in each column, linear interpolation in r; folder film = mean(0.5 - r) over at least %d columns; window film = mean over folders minus the wall shift." % (CLEAR, MIN_COLS),
     "- smear = 100 sum((1-alpha) V, alpha > 0.5) / sum((1-alpha) V).",
     "- body region = 4-connected region (scipy.ndimage.label, cross structure) of cells with alpha < 0.5 containing the nose cell (the highest-x axis-row cell with alpha < 0.5). tail_b = smallest cell-centre x of body cells in the axis row; L_b = nose (interpolated, item 1) - tail_b.",
     "- out_b / wake_b = 100 sum((1-alpha) V) over cells NOT in the body region (including all cells with alpha >= 0.5) with x < tail_b - 0.05 (wake_b: and x > tail_b - 2.0), divided by the total gas of the folder.", "",
     "## T3_rhoST per folder (%d folders)" % len(rs), ""] + table(rs) + [""]
L += ["## T3_prod per window folder (%d folders)" % len(rp), ""] + table(rp) + [""]


def srow(name, key, fmt=f6):
    return "| %s | %s | %s |" % (name, fmt(ss[key]), fmt(sp[key]))


L += ["## Summary", "", "| quantity | T3_rhoST | T3_prod |", "|---|---|---|",
      "| folders in window | %d | %d |" % (ss["nfold"], sp["nfold"]),
      "| n used for Fr, out_b, wake_b, L_b (folders with a nose) | %d | %d |" % (ss["n"], sp["n"]),
      "| folders without any axis cell with alpha < 0.5 (nose undefined) | %s | %s |" % (", ".join(ss["no_nose"]) or "none", ", ".join(sp["no_nose"]) or "none"),
      "| Fr (slope of nose vs t) | %s | %s |" % (f6(ss["Fr"]), f6(sp["Fr"])),
      "| Fr standard error | %s | %s |" % (f6(ss["se"]), f6(sp["se"])),
      "| film, mean of 0.5 - r over folders (before wall shift) | %s | %s |" % (f6(ss["film_raw"]), f6(sp["film_raw"])),
      "| film (minus wall shift) | %s | %s |" % (f6(ss["film"]), f6(sp["film"])),
      "| folders with a valid film (>= 3 columns) | %d | %d |" % (ss["nfilm"], sp["nfilm"]),
      srow("smear mean (%)", "smear"), srow("out_b mean (%)", "out_b"), srow("wake_b mean (%)", "wake_b"),
      srow("L_b mean", "Lb_mean"), srow("L_b min", "Lb_min"), srow("L_b max", "Lb_max"),
      "| L_b min at t / max at t | %s / %s | %s / %s |" % (f6(ss["Lb_min_t"]), f6(ss["Lb_max_t"]), f6(sp["Lb_min_t"]), f6(sp["Lb_max_t"])), "",
      "- shift = 100 (Fr_rhoST / Fr_prod - 1) = %s %%" % f6(shift),
      "- drift of total gas in T3_rhoST = 100 (G_last / G_first - 1) = %s %% (G at t = %s: %s; G at t = %s: %s)" % (f6(drift), g_win[0][1], f6(G_first), g_win[-1][1], f6(G_last)), "",
      "## Total gas G = sum((1-alpha) V) of T3_rhoST at every folder", "", "| t (s) | G | change vs previous folder (%) |", "|---|---|---|"]
for i, (t, nm, G) in enumerate(gas_rows):
    prev = gas_rows[i - 1][2] if i > 0 else None
    L.append("| %s | %s | %s |" % (nm, f6(G), "" if prev is None else f6(100.0 * (G / prev - 1.0))))
L += ["", "Largest single-step change among folders 0.2 to 8.4: %s %%; minimum %s %%." % (f6(np.abs(dG).max()), f6(dG.min())), "",
      "## Diagnostics (body region)", "", "### T3_rhoST", ""] + diag_table(rs) + ["", "### T3_prod", ""] + diag_table(rp) + [""]
# nose monotonic checks
for nm_, rows_ in (("T3_rhoST", rs), ("T3_prod", rp)):
    okr = [r for r in rows_ if np.isfinite(r["nose"])]
    dn = np.diff([r["nose"] for r in okr]); dl = np.diff([r["L_b"] for r in okr])
    L.append("- %s: smallest nose step between consecutive folders with a nose %s, largest |L_b step| %s." % (nm_, f6(dn.min()), f6(np.abs(dl).max())))
(HERE / "checks/B9_t3check.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))

# --- comparison with the given T3_prod values, added only after the file above was written ---
cmp_ = [("Fr", sp["Fr"], EXPECT["Fr"], 1e-5), ("Fr standard error", sp["se"], EXPECT["Fr_se"], 1e-5), ("n", float(sp["n"]), float(EXPECT["n"]), 0.0),
        ("film", sp["film"], EXPECT["film"], 1e-4), ("out_b (%)", sp["out_b"], EXPECT["out_b"], 0.01),
        ("wake_b (%)", sp["wake_b"], EXPECT["wake_b"], 0.01), ("smear (%)", sp["smear"], EXPECT["smear"], 0.01)]
C = ["", "## Comparison of T3_prod with the given values", "", "| quantity | here | given | difference | tolerance | result |", "|---|---|---|---|---|---|"]
allok = True
for name, got, want, tol in cmp_:
    ok = bool(abs(got - want) <= tol); allok &= ok
    C.append("| %s | %s | %s | %+.3e | %g | %s |" % (name, f6(got), f6(want), got - want, tol, "PASS" if ok else "FAIL"))
C += ["", "Tolerances as given: Fr 1e-5, film 1e-4, shares 0.01 percentage points. No expected L_b was given, so L_b is not compared. The standard error and n are compared with tolerance 1e-5 and exact, respectively (my choice, not given).", "",
      "## Result", "", "**%s**" % ("PASS" if allok else "FAIL")]
with open(HERE / "checks/B9_t3check.md", "a", encoding="utf-8") as fh:
    fh.write("\n".join(C) + "\n")
print("\n".join(C))
