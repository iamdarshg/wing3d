import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f5
from wing3d.solver import solve, assemble

mesh, wakes, sref, info = build_f5(scale=0.1)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(2.0)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
C = mesh.centroid
for sgn, tag in [(1, 'up'), (-1, 'lo')]:
    m = (np.abs(C[:, 1] - 0.30) < 0.06) & ((C[:, 2] * sgn) > 0.0005)
    print(tag, 'n=', m.sum())
    if m.sum():
        o = np.argsort(C[m, 0])
        x, cp = C[m, 0][o], res['cp'][m][o]
        print('  x range:', x.min().round(3), x.max().round(3))
        print('  Cp[:8]:', np.round(cp[:8], 2))
