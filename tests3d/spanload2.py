import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=16)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
dy = 6.0 / ns
print('npl=', npl, 'total body panels should be', ns * npl)
for j in range(ns):
    ids = [j * npl + i for i in range(npl)]
    Fz = (-res['cp'][ids, None] * mesh.normal[ids]
          * mesh.area[ids, None]).sum(axis=0)[2]
    ymid = float(np.mean(mesh.centroid[ids][:, 1]))
    print('strip', j, 'y=%.2f' % ymid, 'sec_cl=', round(float(Fz / dy), 4))
