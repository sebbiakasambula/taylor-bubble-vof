"""Run interFoam on one case in Docker, in parallel when the mesh is big enough to need it.

  python run_case.py <case folder>            1 process up to 20 000 cells, 2 above that
  python run_case.py <case folder> --np 4     force the number of processes

So here the rule is simple. Up to 20 000 cells a case runs on 1 process, above that on 2, one per physical core. The
laptop has 2 physical cores and 4 threads, and CFD gets little from hyperthreads, so 2 is the efficient number, and it
leaves the other 2 threads free so the computer stays usable during long runs. The cut is along the pipe axis only (decomposition (N 1 1)), since the
mesh is a one-cell wedge, and the pieces are put back together with reconstructPar at the end, so the analysis reads
the same folders as for a serial run. The cell count comes from log.checkMesh, so run the case setup first.
"""
import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

IMAGE = "opencfd/openfoam-run:2412"
CELLS_PER_PROCESS = 20000
HDR = "FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\n"


def usable_threads():
    return max((os.cpu_count() or 2) // 2, 1)        # physical cores


def cell_count(case):
    return int((Path(case) / "log.checkMesh").read_text().split("cells:")[1].split()[0])


def processes(cells, limit=None):
    """How many processes a mesh of this size gets, one per 20 000 cells, up to the physical cores."""
    return min(math.ceil(cells / CELLS_PER_PROCESS), limit or usable_threads())


def decompose_dict(n):
    return HDR + f"""
numberOfSubdomains {n};
method          simple;
coeffs
{{
    n           ({n} 1 1);
}}
"""


def foam(case, cmd):
    inner = f"source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1; {cmd}"
    return subprocess.run(["docker", "run", "--rm", "-u", "1000:1000", "-v", f"{Path(case).resolve()}:/data",
                           "-w", "/data", "--entrypoint", "bash", IMAGE, "-c", inner]).returncode


def run(case, np=None):
    case = Path(case)
    n = np or processes(cell_count(case))
    if n == 1:
        return foam(case, "interFoam > log.interFoam 2>&1")
    (case / "system" / "decomposeParDict").write_text(decompose_dict(n))
    steps = ("decomposePar -force > log.decomposePar 2>&1",
             f"mpirun --oversubscribe -np {n} interFoam -parallel > log.interFoam 2>&1",
             "reconstructPar > log.reconstructPar 2>&1")
    for cmd in steps:
        code = foam(case, cmd)
        if code != 0:
            print(f"failed: {cmd} (exit {code})")
            return code
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--np", type=int)
    a = ap.parse_args()
    sys.exit(run(a.case, a.np))
