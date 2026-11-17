import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import Normalize; import matplotlib.cm as cm
from config import PV, WORK
R=str(WORK)+'/'; P=str(PV)+'/strips/'; OUT=R+'fig_t3_fields.pdf'   # strips from paraview/render_flat.py
plt.rcParams.update({'font.size':8.5})
fig=plt.figure(figsize=(6.3,2.75))
g=fig.add_gridspec(2,2,width_ratios=[1,0.022],hspace=0.42,wspace=0.03)
for k,(n,t) in enumerate([('T3_OF','(a) T3, OpenFOAM, $t$ = 20.3 s'),('T3_F','(b) T3, Fluent, $t$ = 20.3 s')]):
    a=fig.add_subplot(g[k,0]); a.imshow(mpimg.imread(P+n+'.png'),extent=(-5.325,0.725,-0.55,0.55),aspect='equal'); a.set_title(t,fontsize=8.5,loc='left'); a.set_ylabel('$r/D$')
    if k==1: a.set_xlabel('Distance from the nose, $x/D$')
    else: a.set_xticklabels([])
c=fig.add_subplot(g[:,1]); fig.colorbar(cm.ScalarMappable(Normalize(0,2),'viridis'),cax=c,extend='max',ticks=[0,0.5,1,1.5,2]).set_label('Speed (m/s)')
fig.savefig(OUT,bbox_inches='tight',dpi=300); fig.savefig(R+'fig_t3_fields_preview.png',bbox_inches='tight',dpi=110)
