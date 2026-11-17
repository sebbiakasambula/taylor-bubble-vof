"""Per-snapshot metrics for a Taylor bubble run (OpenFOAM or Fluent cell-centre data).

Pure functions on cell-centre arrays. Axes: x along the pipe (rise direction), y = radius, lengths in D.
alpha is the liquid fraction, so the gas fraction is g = 1 - alpha.

Two tails are used. tail_ax is the AXIS tail (tb_post.tail_position): walk down the axis from the nose to the
first liquid. It is the rule behind the verified film numbers. x_low is the BODY tail: the lowest cell-centre x
of the 4-connected region of alpha < 0.5 cells that contains the nose cell. It needs a structured (x, r) grid.
"""
import math
import sys
import numpy as np
from scipy import ndimage

# tb_post.py is in this folder
import tb_post

KEYS = ("nose", "tail_ax", "film_raw", "ncol", "film", "L_ax", "G", "out_ax", "wake_ax", "smear",
        "pk_all", "pk_liq", "pk_gas", "pk_x", "pk_r", "pk_alpha", "gap05",
        "x_low", "L_body", "skirt", "main", "out_body", "wake_body", "pk_region",
        "u_core", "fast_core", "wk_nfast_ax", "wk_fastgas_ax", "wk_nfast_body", "wk_fastgas_body",
        "tail_b", "L_b", "skirt_b", "out_b", "wake_b", "wk_nfast_b", "wk_fastgas_b")

WALL_OPENFOAM = 0.5 * math.cos(math.radians(2.5))  # wedge mid-plane wall, 0.499524111
WALL_FLUENT = 0.5


def structured_index(x, y):
    """Grid indices (ix, ir) of every cell, from the unique x and y rounded to 6 decimals. Raises if the cells
    do not fill a structured x-r grid without collisions."""
    ux, ix = np.unique(np.round(x, 6), return_inverse=True)
    uy, ir = np.unique(np.round(y, 6), return_inverse=True)
    ix, ir = ix.ravel(), ir.ravel()
    if len(ux) * len(uy) != len(x) or len(np.unique(ix * len(uy) + ir)) != len(x):
        raise ValueError("cells do not form a structured grid: %d x-values, %d r-values, %d cells" % (len(ux), len(uy), len(x)))
    return len(ux), len(uy), ix, ir


def body_region(x, y, alpha):
    """Boolean per cell: the 4-connected region of alpha < 0.5 that contains the nose cell (the highest-x cell
    on the axis row with alpha < 0.5). None if there is no gas on the axis."""
    nx, ny, ix, ir = structured_index(np.asarray(x, float), np.asarray(y, float))
    gas = np.zeros((nx, ny), bool)
    gas[ix, ir] = np.asarray(alpha, float) < 0.5
    axis_gas = np.flatnonzero(gas[:, 0])
    if len(axis_gas) == 0:
        return None
    lab, _ = ndimage.label(gas)  # default structure: 4-connectivity
    inose = axis_gas[-1]
    return (lab == lab[inose, 0])[ix, ir]


def _wake_fast(x, alpha, U_mag, edge):
    """Cells in the 2 D below `edge` (edge - 2.0 < x < edge - 0.05) faster than 0.3 m/s, and the percentage of them
    that hold gas (alpha <= 0.999)."""
    m = (x > edge - 2.0) & (x < edge - 0.05) & (U_mag > 0.3)
    n = int(m.sum())
    return n, (100.0 * float((m & (alpha <= 0.999)).sum()) / n if n else float("nan"))


def snapshot_metrics(x, y, alpha, U_mag, V, wall, clear):
    """Metrics of one snapshot. NaN where a quantity is not defined (no nose, no tail, fewer than 3 columns)."""
    nan = float("nan")
    x = np.asarray(x, float)
    y = np.round(np.asarray(y, float), 9)
    alpha = np.asarray(alpha, float)
    U_mag = np.asarray(U_mag, float)
    V = np.asarray(V, float)

    out = {k: nan for k in KEYS}
    out["ncol"] = 0
    out["pk_region"] = "none"

    g = 1.0 - alpha
    gV = g * V
    G = float(gV.sum())
    out["G"] = G

    # peak speeds
    ipk = int(np.argmax(U_mag))
    out["pk_all"] = float(U_mag[ipk])
    liq = alpha > 0.99
    gas = alpha < 0.01
    out["pk_liq"] = float(U_mag[liq].max()) if liq.any() else nan
    out["pk_gas"] = float(U_mag[gas].max()) if gas.any() else nan
    out["pk_r"] = float(y[ipk])
    out["pk_alpha"] = float(alpha[ipk])

    # smeared gas: gas held in cells that are less than half gas
    out["smear"] = 100.0 * float(gV[g < 0.5].sum()) / G if G > 0 else nan

    nose = tb_post.nose_position(x, y, alpha)
    tail = tb_post.tail_position(x, y, alpha)
    if nose is not None:
        out["nose"] = float(nose)
        out["pk_x"] = float(x[ipk] - nose)
    if tail is not None:
        out["tail_ax"] = float(tail)
    if nose is None:
        return out

    xi, ri = tb_post.interface_radius(x, y, alpha)

    # radial gap between the interface and the wall half a diameter behind the nose
    if len(xi):
        k = int(np.argmin(np.abs(xi - (nose - 0.5))))
        out["gap05"] = float(wall - ri[k])

    if tail is not None:
        out["L_ax"] = float(nose - tail)

        # film over the developed region, clear of the tail and the nose
        w = (xi > tail + 0.5) & (xi < nose - clear)
        out["ncol"] = int(w.sum())
        if out["ncol"] >= 3:
            out["film_raw"] = float(np.mean(0.5 - ri[w]))
            out["film"] = out["film_raw"] - (0.5 - wall)

        # gas shares below the axis tail, in percent of the gas in the whole domain
        if G > 0:
            out["out_ax"] = 100.0 * float(gV[x < tail - 0.05].sum()) / G
            out["wake_ax"] = 100.0 * float(gV[(x > tail - 2.0) & (x < tail - 0.05)].sum()) / G
        out["wk_nfast_ax"], out["wk_fastgas_ax"] = _wake_fast(x, alpha, U_mag, tail)

    # body tail and the quantities that use it
    body = body_region(x, y, alpha)
    if body is not None:
        x_low = float(x[body].min())
        out["x_low"] = x_low
        out["L_body"] = float(nose - x_low)
        if tail is not None:
            out["skirt"] = float(tail - x_low)
        if G > 0:
            out["main"] = 100.0 * float(gV[body].sum()) / G
            out["out_body"] = 100.0 * float(gV[x < x_low - 0.05].sum()) / G
            out["wake_body"] = 100.0 * float(gV[(x > x_low - 2.0) & (x < x_low - 0.05)].sum()) / G
        if body[ipk]:
            out["pk_region"] = "bubble"
        elif x[ipk] < x_low:
            out["pk_region"] = "wake"
        elif x[ipk] > nose:
            out["pk_region"] = "ahead"
        else:
            out["pk_region"] = "beside"
        # gas core: gas-volume-weighted speed over the body cells that are pure gas
        core = body & (alpha < 0.01)
        gc = gV[core]
        if gc.sum() > 0:
            uc = U_mag[core]
            out["u_core"] = float((gc * uc).sum() / gc.sum())
            out["fast_core"] = 100.0 * float(gc[uc > 2.0].sum() / gc.sum())
        out["wk_nfast_body"], out["wk_fastgas_body"] = _wake_fast(x, alpha, U_mag, x_low)

        # tail_b: lowest body cell in the axis row. Steps over isolated alpha >= 0.5 axis cells inside the
        # bubble (the body connects around them) and ignores wall slivers. Detached gas = gas outside the body.
        tail_b = float(x[body & (y <= y.min() + 1e-9)].min())
        out["tail_b"] = tail_b
        out["L_b"] = float(nose - tail_b)
        out["skirt_b"] = float(tail_b - x_low)
        if G > 0:
            det = ~body
            out["out_b"] = 100.0 * float(gV[det & (x < tail_b - 0.05)].sum()) / G
            out["wake_b"] = 100.0 * float(gV[det & (x > tail_b - 2.0) & (x < tail_b - 0.05)].sum()) / G
        out["wk_nfast_b"], out["wk_fastgas_b"] = _wake_fast(x, alpha, U_mag, tail_b)
    return out
