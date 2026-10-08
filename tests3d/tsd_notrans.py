import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.spatial import cKDTree
from scipy.optimize import newton_krylov
from tsd_kroot import build_case

s, base_jump, base_G = build_case(0.3, 1.25, True)
s.wing_jump = base_jump
s.wake_scale = 1.0
# kill transpiration (test if it pollutes TE Kutta)
s._apply_wing_bc = lambda R: R
nphi = s.nx * s.ny * s.nz


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


sol = newton_krylov(func, s.phi.ravel().copy(), iter=60, verbose=False,
                    f_tol=1e-4, f_rtol=1e-10)
s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
L = s.loads()
jmid = int(np.argmin(np.abs(s.yc)))
ite = s.wing_ix[-1]
dcp = L['cpl'][ite, jmid] - L['cpu'][ite, jmid]
print('no-transpiration linear: CL=%.4f dcp_te=%+.5f (panel 0.115)' % (
    L['CL'], dcp))
