import sys
from pathlib import Path
from config import RUNS
from tb_case import write_static_case
from run_case import foam
c=RUNS/"D_coarseAR1"
info=write_static_case(c,Eo=10.,Nf=2000.,rho_r=1000.,nx=96,nr=12,grading=1.0,n_capillary=1.0/0.3955,n_snap=5)
for cmd in ("blockMesh > log.blockMesh 2>&1","checkMesh > log.checkMesh 2>&1","setFields > log.setFields 2>&1"):
    assert foam(c,cmd)==0,cmd
print("max_dt",info["max_dt"],flush=True)
print("exit",foam(c,"interFoam > log.interFoam 2>&1"),flush=True)
