import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for tg in [0.004, 0.002, 0.001, 0.0005]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12,
                      te_gap=tg)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print('te_gap=%.4f CL=%.4f CDp=%.5f' % (tg, f['CL'], f['CDp']))
