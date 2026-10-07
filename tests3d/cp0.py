import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from wing3d.geometry import build_wing, naca4_section
from wing3d.solver import solve, build_wake_from_meta, assemble

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=32, n_span=10)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
res = solve(mesh, wakes, [1, 0, 0], A=A, Brow=Brow, kutta_mode='doublet')
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
j = ns // 2
print('3D midspan upper Cp(x) at alpha=0:')
for i in range(0, 16, 3):
    print('  x=%.3f cp=%+.4f' % (mesh.centroid[j * npl + i][0],
                                 res['cp'][j * npl + i]))
