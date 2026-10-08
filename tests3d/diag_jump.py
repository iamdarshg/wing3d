import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble,
                           _surface_gradient)
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
# recompute Cp with scaled jump correction
for k in [1.0, 0.5, 0.0]:
    vel = res['vel'] - 0.5 * _surface_gradient(mesh, res['mu'])
    # res['vel'] already includes +0.5 grad; remove then re-add k*
    vel = vel + k * 0.5 * _surface_gradient(mesh, res['mu'])
    vmag = np.linalg.norm(V)
    cp = 1.0 - np.sum(vel ** 2, axis=1) / vmag ** 2
    # re-apply TE averaging
    for w in wakes:
        iu, il = int(w.upper_te), int(w.lower_te)
        m = 0.5 * (cp[iu] + cp[il])
        cp[iu] = m
        cp[il] = m
    f = pressure_forces(mesh, cp, V, sref=6.0)
    print('jump_corr=%.1f CL=%.4f CDp=%.5f' % (k, f['CL'], f['CDp']))
