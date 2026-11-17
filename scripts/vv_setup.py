"""Write the OpenFOAM verification cases and check their meshes. No solver is run.

  V1a (= Fluent S)  static bubble, gravity off, at the worst spurious-current corner
  V3                time step at T1: the capillary cap halved (dt_factor 0.5), medium mesh
  V4                grid family at T1, refinement ratio 1.4 in both directions
  V5                domain length at T1: pipe 30 D instead of 20 D, same cell size

All at T1 (Eo 300, Nf 1000) since the re-plan of 27 Oct 2026; V1b and the Nf 2000 hard point are dropped.
Each case is meshed with blockMesh, checked with checkMesh, and initialised with setFields, which takes
seconds.
"""
import json
import subprocess
import sys
from pathlib import Path

from mesh_design import cells_across, radial_cells
from tb_case import write_case, write_static_case
from tb_references import film_llewellin

from config import RUNS
IMAGE = "opencfd/openfoam-run:2412"
GRID = {"coarse": (69, 343), "medium": (96, 480), "fine": (134, 672)}   # (radial, axial), ratio 1.4
T1 = dict(Eo=300.0, Nf=1000.0, rho_r=1000.0)


def foam(case, cmd):
    inner = f"source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1; {cmd} > log.{cmd.split()[0]} 2>&1"
    return subprocess.run(["docker", "run", "--rm", "-u", "1000:1000", "-v", f"{case}:/data", "-w", "/data",
                           "--entrypoint", "bash", IMAGE, "-c", inner]).returncode


def check(case):
    for cmd in ("blockMesh", "checkMesh", "setFields"):
        if foam(case, cmd) != 0:
            return dict(ok=False, failed=cmd)
    log = (case / "log.checkMesh").read_text()
    cells = int(log.split("cells:")[1].split()[0])
    ok = "Mesh OK" in log
    sel = [int(l.split("Selected")[1].split("/")[0]) for l in (case / "log.setFields").read_text().splitlines()
           if "Selected" in l]
    return dict(ok=ok, cells=cells, gas_cells_selected=sel)


results = {}
print("V1, static bubble")
for tag, Eo, Nf in (("V1a", 10.0, 2000.0),):
    c = RUNS / f"{tag}_static"
    info = write_static_case(c, Eo, Nf, 1000.0)
    r = check(c)
    results[tag] = dict(info={k: v for k, v in info.items() if k != "fluid"}, **r)
    print(f"  {tag} Eo={Eo:g} Nf={Nf:g}: mesh OK={r['ok']}, {r.get('cells')} cells, "
          f"setFields selected {r.get('gas_cells_selected')}; Laplace jump {info['laplace_jump']:.4f}, "
          f"capillary time {info['t_capillary']:.3f}, run to t = {info['end_time']:.2f}, "
          f"capillary velocity {info['capillary_velocity']:.1f}")

d = 0.5 * film_llewellin(T1["Nf"])                       # Llewellin eq. (4.2) gives film over the radius
print(f"\nV3 to V5 at T1, Eo={T1['Eo']:g}, Nf={T1['Nf']:g} (Llewellin film {d:.4f} D)")
print(f"  {'case':<10} {'radial':>6} {'axial':>6} {'cells':>7} {'in film':>8} {'max dt':>9}  checks")
cases = {f"V4_{name}": dict(nr=nr, nx=nx) for name, (nr, nx) in GRID.items()}
cases["V3_halfdt"] = dict(nr=96, nx=480, dt_factor=0.5)
cases["V5_long"] = dict(nr=96, nx=720, L=30.0)            # same axial cell, 50% longer pipe
for tag, kw in cases.items():
    c = RUNS / tag
    info = write_case(c, **kw, **T1)
    r = check(c)
    n_film = cells_across(radial_cells(kw["nr"], 0.25)[0], d)
    results[tag] = dict(kw, cells_in_film=n_film, max_dt=info["max_dt"], end_time=info["end_time"], **r)
    print(f"  {tag:<10} {kw['nr']:6d} {kw['nx']:6d} {r.get('cells', 0):7d} {n_film:8d} {info['max_dt']:9.2e}  "
          f"mesh OK={r['ok']}, setFields {r.get('gas_cells_selected')}")

(RUNS / "vv_setup.json").write_text(json.dumps(results, indent=1))
sys.exit(0 if all(v.get("ok") for v in results.values()) else 1)
