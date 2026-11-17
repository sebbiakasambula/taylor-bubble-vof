import sys,re
from pathlib import Path
import numpy as np
from config import RUNS, WORK
exec(open(Path(__file__).resolve().with_name("snap.py")).read().split("case=Path")[0])
for name in sys.argv[1:]:
    case=RUNS/name
    C=vec(case/'0/C'); V=scal(case/'0/V'); x,y=C[:,0],np.round(C[:,1],9)
    ts=sorted([d.name for d in case.iterdir() if re.fullmatch(r'[0-9.]+',d.name) and d.name!='0'],key=float)
    a0=scal(case/'0/alpha.water'); g0=tb_post.gas_volume(a0,V)
    R=[]
    for t in ts:
        a=scal(case/t/'alpha.water'); U=vec(case/t/'U'); sp=np.linalg.norm(U,axis=1)
        n=tb_post.nose_position(x,y,a); tl=tb_post.tail_position(x,y,a); cx=tb_post.centroid_x(x,a,V)
        xi,ri=tb_post.interface_radius(x,y,a); w=(xi>tl+0.5)&(xi<n-1.5) if tl else np.zeros(len(xi),bool)
        fl=(0.5-ri[w]).mean() if w.sum()>=3 else np.nan
        gas=a<0.01; liq=a>0.99
        # gas below the tail (shed)
        shed=((1-a)*V)[(x<(tl or 0)-0.05)].sum()/((1-a)*V).sum() if tl else np.nan
        R.append((float(t),n,tl,cx,fl,sp.max(),sp[liq].max(),sp[gas].max() if gas.any() else 0,tb_post.gas_volume(a,V),shed))
    R=np.array(R); np.save(WORK/(name+'.npy'),R)
    print(f"== {name} cells {len(x)} snapshots {len(R)} gas drift max {100*np.max(np.abs(R[:,8]/g0-1)):.4f}%")
    for r in R[::4]: print("  t %6.2f nose %6.3f tail %6.3f cen %6.3f film %.4f umax %.2f liq %.2f gas %.2f shed %.3f"%(r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[9]))
    T=R[:,0]
    for lo in np.arange(0,T[-1],5):
        k=(T>=lo)&(T<lo+5)
        if k.sum()>=3: print("   %4.0f-%4.0f s Fr nose %.4f centroid %.4f"%(lo,lo+5,np.polyfit(T[k],R[k,1],1)[0],np.polyfit(T[k],R[k,3],1)[0]))
