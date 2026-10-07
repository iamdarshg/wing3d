import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=16)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
res = solve(mesh, wakes, [np.cos(a), 0, np.sin(a)], A=A, Brow=Brow)
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
ys, cls, muws = [], [], []
dy = 6.0 / ns
for j in range(ns):
    ids = [j * npl + i for i in range(npl)]
    F = (-res['cp'][ids, None] * mesh.normal[ids]
         * mesh.area[ids, None]).sum(axis=0)
    ys.append(mesh.centroid[ids][:, 1].mean())
    cls.append(float(2 * F[2] / dy))
    muws.append(float(res['muw'][j]))
ys = np.array(ys)
plt.figure()
plt.plot(ys, cls, 'o-', label='body sectional cl')
plt.plot(ys, np.array(muws) * 2, 's-', label='2*muw (circ equiv)')
plt.xlabel('span y')
plt.ylabel('cl')
plt.legend()
plt.savefig('loading.png')
F = (-res['cp'][:, None] * mesh.normal * mesh.area[:, None]).sum(axis=0)
print('CL=', round(float(F[2] / 6), 4), 'CDp=', round(float(F[0] / 6), 5))
print('saved loading.png')
