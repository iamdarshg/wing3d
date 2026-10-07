"""Calibrate wing_jump sign: panel-mu-exact bound sheet, linear TSD."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.spatial import cKDTree

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(1.25)
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
C = mesh.centroid
mu = np.asarray(res['mu'])
up_idx = C[:, 2] > 0.001
lo_idx = C[:, 2] < -0.001
Tu, Tl = cKDTree(C[up_idx][:, :2]), cKDTree(C[lo_idx][:, :2])
mu_up, mu_lo = mu[up_idx], mu[lo_idx]


def build_jump(s, sign):
    xc = s.xc[s.wing_ix]
    yc = s.yc[s.wing_jy]
    XX, YY = np.meshgrid(xc, yc, indexing='ij')
    pts = np.stack([XX.ravel(), YY.ravel()], axis=-1)
    iu = Tu.query(pts)[1]
    il = Tl.query(pts)[1]
    jump = (mu_lo[il] - mu_up[iu]).reshape(XX.shape)
    return sign * jump


for sign in (+1.0, -1.0):
    s = TSDSolver(Minf=0.3, alpha_deg=1.25, nx=80, ny=32, nz=32,
                  omega=1.5, wake_gamma=g, linear=True)
    s.wing_jump = build_jump(s, sign)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    s.wake_scale = 1.0
    s.solve(itmax=250, tol=0, verbose=False, ramp=0, init_linear=False)
    L = s.loads()
    print('sign=%+.0f: CL=%.4f (panel descaled ~0.115)' % (sign, L['CL']),
          flush=True)
