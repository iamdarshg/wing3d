import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for nc, ks in [(24, 1.0), (24, -1.0), (48, -1.0)]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=nc, n_span=8)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet', kutta_sign=ks)
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print('nc=', nc, 'ks=', ks, 'CL=%+.4f' % f['CL'],
          'Cpmin=%+.2f' % res['cp'].min())
