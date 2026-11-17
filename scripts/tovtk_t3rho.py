"""T3 snapshots in pipe coordinates (no nose shift): VTK files for ParaView, for the density-weighted-surface-tension run and the production run."""
import numpy as np
from pathlib import Path
exec(open(Path(__file__).with_name('tovtk.py')).read().split("of('V1a_static'")[0])
_write = write
def write(name, x, r, sp, al, shift=0.0):
    # alpha values below the float32 range (e.g. 1e-40) make the ParaView ASCII reader fail and drop the whole file, so flush them to zero
    _write(name, x, r, np.where(np.abs(sp) < 1e-30, 0.0, sp), np.where(np.abs(al) < 1e-30, 0.0, al), shift)
SNAP = [('T3_rhoST', '1', 'T3R_1'), ('T3_rhoST', '2.4', 'T3R_24'), ('T3_rhoST', '5', 'T3R_5'),
        ('T3_rhoST', '8.4', 'T3R_84'), ('T3_prod', '8.272418', 'T3P_827')]
for case, t, name in SNAP:
    of(case, t, name, nose=False)
