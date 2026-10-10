"""E9: Minf-continuation M0.9 -> 0.95 -> 1.0 (branch selection test)."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_coupled import build_case
from scipy.optimize import newton_krylov

s, res, mesh = build_case(0.9, 1.25, False, True)
s.k[:] = 1.0


def solve_epochs(s, epochs=8):
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

    for ep in range(epochs):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenRho = s._face_rho(phix)
        sol = newton_krylov(func, s.phi.ravel(), iter=60, verbose=False,
                            f_tol=1e-5, f_rtol=1e-8)
        s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
    s._frozen = s._frozenA = s._frozenRho = None


def section_cl(s):
    L = s.fp_loads()
    jy = np.asarray(s.wing_jy)
    j0 = jy[np.argmin(np.abs(s.yc[jy]))]
    ix = np.asarray(s.wing_ix)
    xc = s.xc[ix] / s.chord
    cl = float(np.trapezoid((L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc))
    return L['CL'], cl


for M in [0.9, 0.95, 1.0]:
    s.Minf = M
    solve_epochs(s, epochs=6)
    CL, cl = section_cl(s)
    print('M=%.2f: full-wing CL=%.4f centerline cl=%.4f' % (M, CL, cl),
          flush=True)
