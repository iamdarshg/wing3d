import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble,
                           wing_strips_from_structured)
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
nl, ns = mesh.meta['n_loop'], mesh.meta['n_span']
strips, pid = wing_strips_from_structured(nl, ns)
npl = nl - 1
for depth in [1, 2, 3]:
    cp = res['cp'].copy()
    for (up, lo, j) in strips:
        # TE pair + `depth-1` neighbors each side -> pair mean
        idx = []
        for d in range(depth):
            iu = j * npl + d
            il = j * npl + (npl - 1 - d)
            idx += [iu, il]
        mu = np.mean([cp[i] for i in idx[0::2]])
        ml = np.mean([cp[i] for i in idx[1::2]])
        m = 0.5 * (mu + ml)
        for i in idx:
            cp[i] = m
    f = pressure_forces(mesh, cp, V, sref=6.0)
    print('blend depth=%d CL=%.4f CDp=%.5f' % (depth, f['CL'], f['CDp']))
