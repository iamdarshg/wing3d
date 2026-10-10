"""E6: farfield sensitivity -- M1.0 k-frozen with 1.5x box."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_coupled import build_case
from scipy.optimize import newton_krylov

box = {'xlim': (-12, 15), 'ylim': (-10.5, 10.5), 'zlim': (-10.5, 10.5)}
s, res, mesh = build_case(1.0, 1.25, False, True, box)
print('box x:', s.x.min().round(1), s.x.max().round(1),
      'y:', s.y.min().round(1), s.y.max().round(1))
nphi = s.nx * s.ny * s.nz
s.k[:] = 1.0


def func(x):
    s.phi[:] = x.reshape(s.nx, s.ny, s.nz)
    R = s.residual()
    R[0, :, :] = 0
    R[-1, :, :] = 0
    R[:, 0, :] = 0
    R[:, -1, :] = 0
    R[:, :, 0] = 0
    R[:, :, -1] = 0
    return R.ravel()


for epoch in range(10):
    s._frozen = s._sup_mask()
    s._frozenA = s._face_coeff()
    dx = s.dx
    phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
    s._frozenRho = s._face_rho(phix)
    sol = newton_krylov(func, s.phi.ravel(), iter=60, verbose=False,
                        f_tol=1e-5, f_rtol=1e-8)
    s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
    L = s.fp_loads()
    print('  epoch %d: CL=%.4f nsup=%d' % (
        epoch, L['CL'], int(s._sup_mask().sum())), flush=True)
s._frozen = s._frozenA = s._frozenRho = None
L = s.fp_loads()
jy = np.asarray(s.wing_jy)
j0 = jy[np.argmin(np.abs(s.yc[jy]))]
ix = np.asarray(s.wing_ix)
xc = s.xc[ix] / s.chord
cl = float(np.trapezoid((L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc))
print('bigbox: full-wing CL=%.4f centerline cl=%.4f (base 0.1129/0.1157)' % (
    L['CL'], cl))
