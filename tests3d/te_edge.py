import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import build_wake_from_meta
from wing3d.panel import biot_savart

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=16, n_span=6)
wakes = build_wake_from_meta(mesh)
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
j = 3
# TE edge endpoints (upper panel edge, lower panel edge, wake leading edge)
up = mesh.verts[np.array(mesh.faces[j * npl + 0])]
lo = mesh.verts[np.array(mesh.faces[j * npl + npl - 1])]
# find TE edge: upper panel edge containing TE point (x~1)
te = np.array(mesh.meta['te_points'][j])
print('upper panel verts x:', up[:, 0].round(3))
print('lower panel verts x:', lo[:, 0].round(3))
w = wakes[j].wake_panels[0]
print('wake panel0 x:', w[:, 0].round(3))
