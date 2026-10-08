import numpy as np
import os
import time
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces
mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
t0 = time.time()
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
tA = time.time() - t0
t0 = time.time()
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
tS = time.time() - t0
f = pressure_forces(mesh, res['cp'], V, sref=6.0)
tag = 'C' if not os.environ.get('WING3D_NO_C') else 'NP'
print('%s: asm=%.1fs solve=%.1fs CL=%.4f CDp=%.5f' % (tag, tA, tS,
                                                     f['CL'], f['CDp']))
np.save('cases/cp_%s.npy' % tag.lower(), res['cp'])
