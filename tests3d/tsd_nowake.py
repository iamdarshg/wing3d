import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble)
from wing3d.tsd import TSDSolver

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
g = {'y': np.array([mesh.centroid[j * 48][1] for j in range(24)]),
     'gamma': np.array(res['muw'])}
s = TSDSolver(Minf=0.3, alpha_deg=4.0, nx=60, ny=24, nz=24, omega=0.8,
              wake_gamma=g)
# NO farfield (phi=0 box). Ramp wake manually.
s.wake_scale = 0.0
for chunk in range(10):
    s.wake_scale = min(1.0, (chunk + 1) / 10)
    s.solve(itmax=50, tol=1e-9, verbose=False, ramp=0)
    R = s.residual()
    r = float(np.abs(R[1:-1, 1:-1, 1:-1]).max())
    print('chunk', chunk + 1, 'scale=%.1f' % s.wake_scale,
          'max|R|=%.3e' % r, 'finite:', np.isfinite(r), flush=True)
    if not np.isfinite(r):
        break
