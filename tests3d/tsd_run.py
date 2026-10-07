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
s = TSDSolver(Minf=0.3, alpha_deg=4.0, nx=80, ny=32, nz=32, omega=0.7,
              wake_gamma=g)
s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
s.wake_scale = 1.0
nx, ny, nz = s.nx, s.ny, s.nz
ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                         indexing='ij')
red = ((ii + jj + kk) % 2 == 0)
dxm = np.append(s.dx, s.dx[-1])[:, None, None]
dym = np.append(s.dy, s.dy[-1])[None, :, None]
dzm = np.append(s.dz, s.dz[-1])[None, None, :]
D = (2 * max(1 - 0.09, 0.05) / dxm ** 2 + 2 / dym ** 2 + 2 / dzm ** 2)
F = s.farfield
for it in range(300):
    s.phi[0, :, :] = F[0, :, :]
    s.phi[-1, :, :] = F[-1, :, :]
    s.phi[:, 0, :] = F[:, 0, :]
    s.phi[:, -1, :] = F[:, -1, :]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    for mask in (red, ~red):
        R = s.residual()
        upd = np.zeros_like(s.phi)
        upd[mask] = +0.7 * R[mask] / D[mask]
        upd[0, :, :] = 0
        upd[-1, :, :] = 0
        upd[:, 0, :] = 0
        upd[:, -1, :] = 0
        upd[:, :, 0] = 0
        upd[:, :, -1] = 0
        s.phi[mask] += upd[mask]
    if it % 50 == 0:
        R = s.residual()
        print(f'it {it}: max|R|='
              f'{np.abs(R[1:-1, 1:-1, 1:-1]).max():.3e}', flush=True)
L = s.loads()
print('TSD M0.3 full-wake: CL=%.4f (panel 0.367)' % L['CL'])
