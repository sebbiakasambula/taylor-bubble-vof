import re,sys,math
from pathlib import Path
import numpy as np
# tb_post.py is in this folder
import tb_post
def scal(f):
    s=Path(f).read_text(); m=re.search(r'internalField\s+nonuniform List<scalar>\s*(\d+)\s*\(',s)
    n=int(m[1]); return np.fromstring(s[m.end():s.index(')',m.end())],sep=' ',count=n)
def vec(f):
    s=Path(f).read_text(); m=re.search(r'internalField\s+nonuniform List<vector>\s*(\d+)\s*\(',s)
    n=int(m[1]); body=s[m.end():]; body=body[:body.index('\n)')].replace('(',' ').replace(')',' ')
    return np.fromstring(body,sep=' ',count=3*n).reshape(n,3)
case=Path(sys.argv[1]); 
C=vec(case/'0/C'); V=scal(case/'0/V')
sig=float(sys.argv[3]) if len(sys.argv)>3 else 0.0999
rows=[]
ts=sorted([d.name for d in case.iterdir() if re.fullmatch(r'[0-9.]+',d.name) and d.name!='0'],key=float)
for t in ts:
    a=scal(case/t/'alpha.water'); U=vec(case/t/'U'); p=scal(case/t/'p'); sp=np.linalg.norm(U,axis=1)
    gas=a<0.01; liq=a>0.99; itf=~gas&~liq
    Reff=tb_post.effective_radius(a,V); j=tb_post.laplace_jump(C[:,0],C[:,1],a,p,2.0,0.25)
    print(f"t {float(t):6.3f}  max gas {sp[gas].max():7.3f}  max interface {sp[itf].max():7.3f}  max liquid {sp[liq].max():7.4f}  gasV {tb_post.gas_volume(a,V):.6e}  Reff {Reff:.4f}  jump {j if j is None else round(j,4)}  target {2*sig/Reff:.4f}")
