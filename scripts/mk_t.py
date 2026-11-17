import sys,subprocess,re
from pathlib import Path
import numpy as np
from config import RUNS
from tb_case import write_case
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])
base=RUNS
IMAGE="opencfd/openfoam-run:2412"
def foam(case,cmd):
    inner=f"source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1; {cmd}"
    return subprocess.run(["docker","run","--rm","-u","1000:1000","-v",f"{case.resolve()}:/data","-w","/data","--entrypoint","bash",IMAGE,"-c",inner]).returncode
for name,Eo,Nf in (('T2_prod',100.0,50.0),('T3_prod',20.0,500.0)):
    c=base/name
    if c.exists(): print('exists',c); continue
    info=write_case(str(c),Eo=Eo,Nf=Nf,rho_r=1000.0)
    print(name,{k:info[k] for k in info if k in('max_dt','end_time','radius','initial_nose')})
    foam(c,"blockMesh > log.blockMesh 2>&1; checkMesh > log.checkMesh 2>&1; cp -r 0 0.orig 2>/dev/null; setFields > log.setFields 2>&1; postProcess -func writeCellCentres -time 0 >/dev/null 2>&1; postProcess -func writeCellVolumes -time 0 >/dev/null 2>&1")
    a=scal(c/'0/alpha.water'); V=scal(c/'0/V'); fr=((1-a)*V).sum()/V.sum()
    print(name,'cells',len(V),'gas frac',fr,'full gas',fr*15.70796,'B range',fr*15.70796*.995,fr*15.70796*1.005)
