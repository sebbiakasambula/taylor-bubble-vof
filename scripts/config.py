"""Paths used by every script. All roots come from environment variables, with defaults, so the scripts can live
and run anywhere.

P2_RUNS      OpenFOAM run folders, one folder per case (V4_medium, VA_square, T2_prod, T3_prod, ...). Default ~/cfdml/runs_p2.
P2_FLUENT    Fluent exports (T1_Case, T1_square, T1_square_w_completion, T2_case, T3_Case). Default ../Fluent_Data, next to the repo.
P2_WORK      intermediate arrays (.npy), per-snapshot csv files, tables, VTK files, renders and check reports. Default ./work in the repo.
P2_PARAVIEW  the pvbatch executable that runs the scripts in scripts/paraview. Default 'pvbatch'."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

RUNS = Path(os.environ.get("P2_RUNS", Path.home() / "cfdml" / "runs_p2"))
FLUENT = Path(os.environ.get("P2_FLUENT", REPO.parent / "Fluent_Data"))
WORK = Path(os.environ.get("P2_WORK", REPO / "work"))
PARAVIEW = os.environ.get("P2_PARAVIEW", "pvbatch")

PV = WORK / "pv"                    # VTK files from tovtk*.py and prep_film_mesh.py, ParaView strips in pv/strips
CHECKS = WORK / "checks"            # check reports written by the figure scripts
FIGS = REPO / "paper" / "figures"   # finished figures that paper.tex includes; the scripts write to WORK, copy by hand

for _d in (WORK, PV, CHECKS):
    _d.mkdir(parents=True, exist_ok=True)
