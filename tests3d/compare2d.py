import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from wing3d.geometry import build_wing, naca4_section
from wing3d.solver import solve, build_wake_from_meta, assemble
from panel2d import panel_2d

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=32, n_span=10)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
j = ns // 2
xs3, cp3u, cp3l = [], [], []
for i in range(npl):
    x = mesh.centroid[j * npl + i][0]
    z = mesh.centroid[j * npl + i][2]
    if z >= 0:
        xs3.append(x)
        cp3u.append(res['cp'][j * npl + i])
    else:
        cp3l.append(res['cp'][j * npl + i])
print('3D upper: min', round(min(cp3u), 3))
print('3D lower: max', round(max(cp3l), 3))
