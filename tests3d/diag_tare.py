import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
A, Brow = assemble(mesh, [], kutta_mode='doublet')
for adeg in [0, 2, 4, 6, 8]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, [], V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print('wakeless a=%d CL=%.4f CDp=%.5f' % (adeg, f['CL'], f['CDp']))
