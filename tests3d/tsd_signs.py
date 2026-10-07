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
base = tsd_farfield_box(
    TSDSolver(Minf=0.3, alpha_deg=4.0, nx=60, ny=24, nz=24,
              wake_gamma=g),
    wakes_long, np.array(res['muw']))
for ws in [+1.0, -1.0]:
    for fs in [+1.0, -1.0, 0.0]:
        s = TSDSolver(Minf=0.3, alpha_deg=4.0, nx=60, ny=24, nz=24,
                      omega=0.8, wake_gamma=g, wake_sign=ws)
        s.farfield = fs * base if fs != 0.0 else None
        s.wake_scale = 1.0
        s.solve(itmax=120, tol=1e-9, verbose=False, ramp=0)
        R = s.residual()
        r = float(np.abs(R[1:-1, 1:-1, 1:-1]).max())
        cl = s.loads()['CL'] if np.isfinite(r) else float('nan')
        print(f'wake {ws:+.0f} far {fs:+.0f}: max|R|={r:.3e} CL={cl:+.4f}',
              flush=True)
