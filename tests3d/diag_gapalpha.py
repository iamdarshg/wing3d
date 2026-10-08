import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12,
                  te_gap=0.0025)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
for adeg in [0, 2, 4, 6, 8]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    cdi = f['CL'] ** 2 / (np.pi * 6.0 * 0.9)
    print('a=%d CL=%.4f CDp=%.5f CDi_est=%.5f ratio=%.2f' % (
        adeg, f['CL'], f['CDp'], cdi, f['CDp'] / max(cdi, 1e-9)))
