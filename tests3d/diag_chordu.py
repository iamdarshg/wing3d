import numpy as np
import sys
import time
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for nc in [16, 24, 32, 48]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=nc, n_span=12,
                      cosine=False)
    wakes = build_wake_from_meta(mesh)
    t0 = time.time()
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    tA = time.time() - t0
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print('n_chord=%d npan=%d asm=%.0fs CL=%.4f CDp=%.5f' % (
        nc, mesh.npanels, tA, f['CL'], f['CDp']))
