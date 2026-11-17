"""B10 prep: cut the front wedge face of two OpenFOAM meshes down to a window behind the T1 nose and write
legacy VTK files for ParaView, plus an npz with the numbers the figure script needs.

Reads only (never writes) <runs>/<case>/constant/polyMesh, <case>/0/C and <case>/19.700352/alpha.water.
Polygons are the real front-wedge-face cells of the mesh, taken from points, faces, owner. Plot coordinates are
x - x_nose and the true radius r = sqrt(y^2 + z^2) of the front-face vertices (wall at r = 0.5).
Cell-centre coordinates and the film numbers are in the wedge mid-plane (y, wall at 0.5 cos 2.5 deg).
"""
import math, re, sys
from pathlib import Path
import numpy as np

from config import RUNS, PV
OUT = PV
TIME = '19.700352'
CASES = {'V4_medium': 11.3199, 'VA_square': 11.5606}
XW = (-2.56, -2.44)          # window in x - x_nose
RW = (0.42, 0.505)           # window in r
MARG_X, MARG_R = 0.06, 0.012  # extra cells kept around the window so edges run off the panel and the contour has neighbours
A = math.radians(2.5)
Y_WALL = 0.5 * math.cos(A)


def body(txt, key):
    """Text of the first '( ... )' list after the line holding the list length."""
    m = re.search(r'\n(\d+)\s*\n\(\s*\n', txt)
    n = int(m[1]); i = m.end()
    j = txt.index('\n)', i)
    return n, txt[i:j]


def read_points(case):
    t = (RUNS / case / 'constant/polyMesh/points').read_text()
    n, b = body(t, 'points')
    return np.fromstring(b.replace('(', ' ').replace(')', ' '), sep=' ', count=3 * n).reshape(n, 3)


def read_owner(case):
    t = (RUNS / case / 'constant/polyMesh/owner').read_text()
    n, b = body(t, 'owner')
    return np.fromstring(b, sep=' ', count=n, dtype=np.int64)


def front_patch(case):
    t = (RUNS / case / 'constant/polyMesh/boundary').read_text()
    m = re.search(r'\bfront\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);', t)
    return int(m[2]), int(m[1])


def read_faces(case, start, n):
    t = (RUNS / case / 'constant/polyMesh/faces').read_text()
    nf, b = body(t, 'faces')
    lines = b.split('\n')
    assert len(lines) == nf, (len(lines), nf)
    out = []
    for ln in lines[start:start + n]:
        k = ln.index('(')
        out.append([int(v) for v in ln[k + 1:-1].split()])
    return out


def read_scalar(f):
    t = Path(f).read_text()
    m = re.search(r'internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(', t)
    return np.fromstring(t[m.end():t.index(')', m.end())], sep=' ', count=int(m[1]))


def read_vector(f):
    t = Path(f).read_text()
    m = re.search(r'internalField\s+nonuniform List<vector>\s*(\d+)\s*\(', t)
    b = t[m.end():t.index('\n)', m.end())]
    return np.fromstring(b.replace('(', ' ').replace(')', ' '), sep=' ', count=3 * int(m[1])).reshape(-1, 3)


def write_unstructured(fn, pts, quads, cell_scalars):
    with open(fn, 'w') as f:
        f.write('# vtk DataFile Version 3.0\nfilm mesh window\nASCII\nDATASET UNSTRUCTURED_GRID\n')
        f.write(f'POINTS {len(pts)} double\n')
        for p in pts: f.write(f'{p[0]:.9g} {p[1]:.9g} 0\n')
        f.write(f'CELLS {len(quads)} {5 * len(quads)}\n')
        for q in quads: f.write('4 ' + ' '.join(str(int(v)) for v in q) + '\n')
        f.write(f'CELL_TYPES {len(quads)}\n' + '9\n' * len(quads))
        f.write(f'CELL_DATA {len(quads)}\n')
        for name, arr in cell_scalars.items():
            f.write(f'SCALARS {name} double 1\nLOOKUP_TABLE default\n')
            for v in arr: f.write(f'{v:.9g}\n')


def write_structured(fn, xc, rc, alpha2d):
    """Cell-centre grid: points at the cell centres, point data alpha, so the 0.5 contour is the one through the centre values."""
    nx, nr = alpha2d.shape
    with open(fn, 'w') as f:
        f.write('# vtk DataFile Version 3.0\ncell centre grid\nASCII\nDATASET STRUCTURED_GRID\n')
        f.write(f'DIMENSIONS {nx} {nr} 1\nPOINTS {nx * nr} double\n')
        for j in range(nr):
            for i in range(nx): f.write(f'{xc[i]:.9g} {rc[j]:.9g} 0\n')
        f.write(f'POINT_DATA {nx * nr}\nSCALARS alpha.water double 1\nLOOKUP_TABLE default\n')
        for j in range(nr):
            for i in range(nx): f.write(f'{alpha2d[i, j]:.9g}\n')


for case, nose in CASES.items():
    P = read_points(case); own = read_owner(case)
    s0, nfr = front_patch(case)
    F = read_faces(case, s0, nfr)
    cells = own[s0:s0 + nfr]                       # owner cell of each front-face
    alpha = read_scalar(RUNS / case / TIME / 'alpha.water')
    C = read_vector(RUNS / case / '0/C')           # OpenFOAM cell centres (written earlier by the solver set-up)
    ncell = len(alpha); assert C.shape[0] == ncell == own.max() + 1
    # keep quads only that lie in the extended window
    sel = []
    for k, f in enumerate(F):
        if len(f) != 4: continue
        v = P[f]; xr = v[:, 0] - nose; rr = np.hypot(v[:, 1], v[:, 2])
        if xr.max() < XW[0] - MARG_X or xr.min() > XW[1] + MARG_X or rr.max() < RW[0] - MARG_R or rr.min() > RW[1] + MARG_R: continue
        sel.append(k)
    sel = np.array(sel)
    ids = sorted({p for k in sel for p in F[k]})
    remap = {p: i for i, p in enumerate(ids)}
    pts = np.array([[P[p, 0] - nose, math.hypot(P[p, 1], P[p, 2])] for p in ids])
    quads = [[remap[p] for p in F[k]] for k in sel]
    cid = cells[sel]
    a_sel = alpha[cid]
    write_unstructured(OUT / f'film_mesh_{case}.vtk', pts, quads, {'alpha.water': a_sel, 'cell_id': cid.astype(float)})
    # cell geometry in the mid-plane for the numbers
    xc = C[cid, 0] - nose; yc = C[cid, 1]
    # per cell radial edges (mid-plane y of its vertices) and axial edges
    y_lo = np.array([min(P[F[k]][:, 1]) for k in sel]); y_hi = np.array([max(P[F[k]][:, 1]) for k in sel])
    x_lo = np.array([min(P[F[k]][:, 0]) for k in sel]) - nose; x_hi = np.array([max(P[F[k]][:, 0]) for k in sel]) - nose
    yc_formula = (2 / 3) * (y_hi ** 3 - y_lo ** 3) / (y_hi ** 2 - y_lo ** 2)
    print(f'{case}: {len(sel)} quads, {len(ids)} points; max |y_c(C) - formula| = {np.abs(yc - yc_formula).max():.2e}; '
          f'max |z_c| = {np.abs(C[cid, 2]).max():.1e}')
    # cell-centre grid (structured), centres from the solver's own cell-centre file
    ux = np.unique(np.round(xc, 7)); ur = np.unique(np.round(yc, 8))
    assert len(ux) * len(ur) == len(sel), (len(ux), len(ur), len(sel))
    ix = np.searchsorted(ux, np.round(xc, 7)); ir = np.searchsorted(ur, np.round(yc, 8))
    a2 = np.full((len(ux), len(ur)), np.nan); a2[ix, ir] = a_sel
    assert not np.isnan(a2).any()
    write_structured(OUT / f'film_ctr_{case}.vtk', ux, ur / math.cos(A), a2)
    np.savez(OUT / f'film_mesh_{case}.npz', nose=nose, ux=ux, ur=ur, a2=a2, x_lo=x_lo, x_hi=x_hi, y_lo=y_lo, y_hi=y_hi,
             xc=xc, yc=yc, ix=ix, ir=ir)
    print('  columns', len(ux), 'rows', len(ur), 'x range', ux.min(), ux.max(), 'r range', ur.min(), ur.max())
