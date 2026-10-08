"""k-scan by tangency error: does least-squares tangency select k=1?"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box, thickness_slope
from scipy.spatial import cKDTree
from scipy.optimize import newton_krylov
from tsd_kroot import build_case, solve_epoch

s, base_jump, base_G = build_case(0.3, 1.25, True)
s.wake_scale = 1.0
jmid = int(np.argmin(np.abs(s.yc)))
k0 = s.k0
iu = np.asarray(s.wing_ix)
xc = np.clip(s.xc[iu] / s.chord, 0, 1)
st = thickness_slope(xc, s.thick)
wu = st - s.alpha
wl = -st - s.alpha


def tang_rms(s):
    zp = (s.phi[iu[:, None], s.wing_jy, k0 + 1]
          - s.phi[iu[:, None], s.wing_jy, k0]) / s.dz[k0]
    zm = (s.phi[iu[:, None], s.wing_jy, k0]
          - s.phi[iu[:, None], s.wing_jy, k0 - 1]) / s.dz[k0 - 1]
    e = np.concatenate([(zp - wu[:, None]).ravel(),
                        (zm - wl[:, None]).ravel()])
    return float(np.sqrt((e ** 2).mean()))


def func_factory(s):
    def func(x):
        s.phi[:] = x.reshape(s.nx, s.ny, s.nz)
        R = s.residual()
        R[0] = 0
        R[-1] = 0
        R[:, 0] = 0
        R[:, -1] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        return R.ravel()
    return func


func = func_factory(s)
for k in [0.5, 0.75, 1.0, 1.25, 1.5]:
    s.wing_jump = k * base_jump
    s.wake_G = k * base_G
    s.phi[1:-1, 1:-1, 1:-1] = 0
    F = s.farfield
    s.phi[0] = F[0]
    s.phi[-1] = F[-1]
    s.phi[:, 0] = F[:, 0]
    s.phi[:, -1] = F[:, -1]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    x0 = s.phi.ravel().copy()
    sol = newton_krylov(func, x0, iter=60, verbose=False, f_tol=1e-4,
                        f_rtol=1e-10)
    s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
    L = s.loads()
    print('k=%.2f CL=%.4f tangRMS=%.4f' % (k, L['CL'], tang_rms(s)),
          flush=True)
