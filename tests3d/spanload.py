import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=16)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
f = pressure_forces(mesh, res['cp'], V, sref=6.0)
print('total CL:', round(f['CL'], 4))
cp = res['cp']
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
dy = 6.0 / ns
s = 0.0
for j in range(ns):
    ids = [j * npl + i for i in range(npl)]
    Fz = (-cp[ids, None] * mesh.normal[ids] * mesh.area[ids, None]).sum(axis=0)[2]
    s += Fz
    if j in [0, 1, 7, 8, 14, 15]:
        ymid = float(np.mean(mesh.centroid[ids][:, 1]))
        print('strip', j, 'y=%.2f' % ymid, 'sec_cl=', round(float(2 * Fz / dy), 4))
Fztot = float((-cp[:, None] * mesh.normal * mesh.area[:, None]).sum(axis=0)[2])
print('sum strips Fz:', round(s, 4), ' total Fz:', round(Fztot, 4))
