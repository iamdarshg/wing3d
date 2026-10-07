import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                            wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
meta = mesh.meta
nl, ns = meta['n_loop'], meta['n_span']
strips, _ = wing_strips_from_structured(nl, ns)
info = [(up, lo, meta['te_points'][j], meta['te_points'][j + 1])
        for (up, lo, j) in strips]
wakes_long = build_wake(mesh, info, length=9.0, n_panels=6)
g = {'y': np.array([mesh.centroid[j * 48][1] for j in range(24)]),
     'gamma': np.array(res['muw'])}
s = TSDSolver(Minf=0.3, alpha_deg=4.0, nx=80, ny=32, nz=32, omega=0.8,
              wake_gamma=g)
s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
s.wake_scale = 0.0
for chunk in range(10):
    s.wake_scale = min(1.0, (chunk + 1) / 10)
    s.solve(itmax=50, tol=1e-9, verbose=False, ramp=0)
    L = s.loads()
    print('after', (chunk + 1) * 50, 'scale=%.1f' % s.wake_scale,
          ': CL=%.4f' % L['CL'], flush=True)
