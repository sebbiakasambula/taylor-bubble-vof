"""Post-processing for the Taylor bubble: rise velocity, nose position, film thickness, and gas volume
for the mass-conservation check. Works on cell-centre arrays in any order, so it does not depend on how
the fields were read. Tested on a synthetic bubble with known answers in test_tb_post.py.

Two measures of position are provided. The GAS CENTROID is the primary one. It is an integral over the
whole bubble, so it moves smoothly and its velocity is the rise velocity once the shape is steady. The
NOSE is a single interface crossing, useful for the shape but noisier, since along the axis it is
located only to a fraction of an axial cell.

The terminal-velocity window is chosen by a stated rule, never by eye. A local velocity is fitted over
LOCAL_SAMPLES consecutive snapshots; the window starts at the earliest snapshot after which that local
velocity stays within TERMINAL_TOL of the velocity fitted over the final FINAL_FRACTION of the record.
A raw derivative is not used, because differencing positions known only to a fraction of a cell
amplifies that error into the velocity.
"""
import numpy as np

TERMINAL_TOL = 0.01
FINAL_FRACTION = 0.3
LOCAL_SAMPLES = 5


def _axis_row(Cx, Cy):
    inner = Cy <= np.min(Cy) + 1e-9
    o = np.argsort(Cx[inner])
    return Cx[inner][o], np.where(inner)[0][o]


def _crossing(x0, x1, a0, a1):
    return x0 + (0.5 - a0) * (x1 - x0) / (a1 - a0)


def nose_position(Cx, Cy, alpha):
    """Upper crossing of alpha = 0.5 along the axis. None if there is no gas."""
    x, idx = _axis_row(Cx, Cy)
    a = alpha[idx]
    gas = np.where(a < 0.5)[0]
    if len(gas) == 0:
        return None
    k = gas[-1]
    if k + 1 >= len(x):
        return float(x[k])
    return float(_crossing(x[k], x[k + 1], a[k], a[k + 1]))


def tail_position(Cx, Cy, alpha):
    """Lower end of the gas column under the nose, along the axis. None if there is no gas or it touches the bottom.

    I walk down from the nose and stop at the first liquid, so gas shed at the tail is not counted as bubble. Taking
    the lowest gas on the axis instead put the tail on shed gas several D below the bubble and spoiled the film window.
    """
    x, idx = _axis_row(Cx, Cy)
    a = alpha[idx]
    gas = np.where(a < 0.5)[0]
    if len(gas) == 0:
        return None
    k = gas[-1]
    while k > 0 and a[k - 1] < 0.5:
        k -= 1
    if k == 0:
        return None
    return float(_crossing(x[k - 1], x[k], a[k - 1], a[k]))


def centroid_x(Cx, alpha, V):
    """Axial position of the gas centroid."""
    g = (1.0 - alpha) * V
    return float(np.sum(Cx * g) / np.sum(g))


def gas_volume(alpha, V):
    return float(np.sum((1.0 - alpha) * V))


def laplace_jump(Cx, Cy, alpha, p, centre_x, radius, margin=0.3):
    """Static bubble (V1): mean pressure inside minus outside, sampled clear of the interface.

    Cells within margin * radius of the interface are left out on both sides, because the pressure
    there is smeared over the cells the interface occupies. Exact value for a sphere: 2 sigma / R.
    """
    rr = np.hypot(Cx - centre_x, Cy)
    inside = (rr < (1.0 - margin) * radius) & (alpha < 0.01)
    outside = (rr > (1.0 + margin) * radius) & (rr < 2.0 * radius) & (alpha > 0.99)
    if inside.sum() == 0 or outside.sum() == 0:
        return None
    return float(np.mean(p[inside]) - np.mean(p[outside]))


def effective_radius(alpha, V, half_angle_deg=2.5):
    """Radius of the sphere with the bubble's actual gas volume, for the Laplace check.

    setFields marks whole cells, so the initial bubble is not exactly the requested sphere. A flat wedge of
    half-angle theta holds the fraction tan(theta) / pi of the full body of revolution, so the full gas
    volume is V_wedge * pi / tan(theta) and R_eff = (3 V_wedge / (4 tan(theta)))^(1/3).
    """
    v = gas_volume(alpha, V)
    return float((3.0 * v / (4.0 * np.tan(np.radians(half_angle_deg)))) ** (1.0 / 3.0))


def max_speed(U):
    """Largest velocity magnitude; for the static bubble this is the spurious current."""
    return float(np.max(np.linalg.norm(np.atleast_2d(U), axis=1)))


def interface_radius(Cx, Cy, alpha):
    """Radius of the alpha = 0.5 surface in every axial column that contains gas."""
    xr = np.round(Cx, 10)
    o = np.lexsort((Cy, xr))
    xs, rs, al = xr[o], Cy[o], alpha[o]
    starts = np.flatnonzero(np.r_[True, xs[1:] != xs[:-1]])
    ends = np.r_[starts[1:], len(xs)]
    out_x, out_r = [], []
    for s, e in zip(starts, ends):
        a, r = al[s:e], rs[s:e]
        cross = np.flatnonzero((a[:-1] < 0.5) & (a[1:] >= 0.5))
        if len(cross):
            k = cross[-1]
            out_x.append(xs[s])
            out_r.append(_crossing(r[k], r[k + 1], a[k], a[k + 1]))
    return np.array(out_x), np.array(out_r)


def film_thickness(Cx, Cy, alpha, r_wall, nose_clear=1.5, tail_clear=0.5):
    """Mean film thickness over the developed region, clear of the nose and the tail.

    Returns (delta, spread, n_columns). The spread is the standard deviation over the window, and a
    large spread means the film is not developed there.
    """
    nose, tail = nose_position(Cx, Cy, alpha), tail_position(Cx, Cy, alpha)
    if nose is None or tail is None:
        return None, None, 0
    x, r = interface_radius(Cx, Cy, alpha)
    m = (x > tail + tail_clear) & (x < nose - nose_clear)
    if m.sum() < 3:
        return None, None, int(m.sum())
    d = r_wall - r[m]
    return float(d.mean()), float(d.std()), int(m.sum())


def _local_velocity(t, x, k=LOCAL_SAMPLES):
    n, h = len(t), k // 2
    u = np.full(n, np.nan)
    for i in range(n):
        lo, hi = max(0, i - h), min(n, i + h + 1)
        if hi - lo >= 3:
            u[i] = np.polyfit(t[lo:hi], x[lo:hi], 1)[0]
    return u


def rise_velocity(t, pos, tol=TERMINAL_TOL, final_fraction=FINAL_FRACTION, t_start=None):
    """Terminal velocity of a position record, with the fit window chosen by the stated rule.

    If t_start is given, the window is taken from it instead, which is how the noisier nose record is
    fitted over the window found from the centroid (see bubble_velocity).

    Always returns the same keys. If terminal velocity was not reached, 'terminal' is False and the
    velocity is the best estimate available, flagged so that it is never used silently.
    """
    t = np.asarray(t, float)
    pos = np.array([np.nan if p is None else p for p in pos], float)
    ok = np.isfinite(pos)
    t, pos = t[ok], pos[ok]
    out = dict(U_b=np.nan, terminal=False, t_start=np.nan, t_end=np.nan, n_window=0,
               fit_rms=np.nan, discarded_fraction=np.nan, reason="")
    if len(t) < 6:
        out["reason"] = "fewer than 6 usable samples"
        return out
    if t_start is not None:
        start = int(np.searchsorted(t, t_start))
    else:
        u = _local_velocity(t, pos)
        n_final = max(3, int(round(final_fraction * len(t))))
        u_final = np.polyfit(t[-n_final:], pos[-n_final:], 1)[0]
        within = np.abs(u - u_final) <= tol * abs(u_final)
        outside = np.flatnonzero(~within)
        start = 0 if len(outside) == 0 else outside[-1] + 1
        if len(t) - start < 3:
            out.update(U_b=float(u_final), t_end=float(t[-1]),
                       reason="velocity still changing at the end of the record")
            return out
    if len(t) - start < 3:
        out["reason"] = "window too short"
        return out
    p = np.polyfit(t[start:], pos[start:], 1)
    resid = pos[start:] - np.polyval(p, t[start:])
    out.update(U_b=float(p[0]), terminal=True, t_start=float(t[start]), t_end=float(t[-1]),
               n_window=int(len(t) - start), fit_rms=float(np.sqrt(np.mean(resid ** 2))),
               discarded_fraction=float(start / len(t)),
               reason="terminal" if t_start is None else "window taken from the centroid")
    return out


def bubble_velocity(t, centroid, nose, steady_tol=0.01):
    """Rise velocity of the bubble, with the terminal window found from the centroid.

    The window is a property of the bubble, not of how its position is measured, so it is found from
    the centroid, an integral measure, and the nose is then fitted over the same window. Along the axis
    the nose is located only to a fraction of an axial cell, and applied to the nose alone the terminal
    rule fails on noise, not physics. If the nose and centroid velocities agree, the bubble's shape has
    stopped changing, which is a check that comes for free.
    """
    rc = rise_velocity(t, centroid)
    rn = rise_velocity(t, nose, t_start=rc["t_start"] if rc["terminal"] else None)
    steady = bool(rc["terminal"] and rn["terminal"] and abs(rn["U_b"] / rc["U_b"] - 1) < steady_tol)
    return dict(U_b=rc["U_b"], terminal=rc["terminal"], shape_steady=steady, centroid=rc, nose=rn)
