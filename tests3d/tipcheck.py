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
print('strip | y | sec_cl | cap_cl')
for j in [0, 1, 3, 7, 8]:
    ids = [j * npl + i for i in range(npl)]
    F = (-res['cp'][ids, None] * mesh.normal[ids]
         * mesh.area[ids, None]).sum(axis=0)
    print(j, round(float(mesh.centroid[ids][:, 1].mean()), 2),
          round(float(2 * F[2] / dy), 4))
# tip caps
for jcap, name in [(0, 'left-tip-cap'), (1, 'right-tip-cap')]:
    base = ns * npl
    ncap = (nl - 1)
    start = base + (0 if jcap == 0 else ncap)
    ids = list(range(start, start + ncap))
    F = (-res['cp'][ids, None] * mesh.normal[ids]
         * mesh.area[ids, None]).sum(axis=0)
    print(name, 'CL-contrib:', round(float(F[2] / 6.0), 4))
F = (-res['cp'][:, None] * mesh.normal * mesh.area[:, None]).sum(axis=0)
print('total CL/6 =', round(float(F[2] / 6.0), 4))
