"""TSD/fullpot at M=1.0 (sonic freestream): coupled k + panel-target rows.

At Minf=1 the TSD linear term vanishes (purely nonlinear) and PG is
singular -- this exercises the fullpot path where it matters most.
Validates vs the TRUE M1.0 OpenFOAM anchor when it lands.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'tests3d')
from tsd_coupled import build_case, panel_aft_target
from scipy.optimize import newton_krylov


def run_m1(alpha_deg=1.25, fullpot=True, span=6.0):
    from wing3d.tsd import TSDSolver
    s, res, mesh = build_case(1.0, alpha_deg, False, fullpot, None,
                              80, 32, 32, span)
    nphi = s.nx * s.ny * s.nz
    s.k[:] = 1.0
    tare = panel_aft_target(mesh, res, s.yc[s.wing_jy])
    print('panel target mean=%.4f' % tare.mean(), flush=True)

    def func(x):
        s.phi[:] = x[:nphi].reshape(s.nx, s.ny, s.nz)
        s.k[:] = x[nphi:]
        R = s.residual()
        R[0, :, :] = 0
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        Rk = (s.kutta_residual(upstream=2) - tare) * 10.0
        return np.concatenate([R.ravel(), Rk])

    for epoch in range(10):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        if fullpot:
            dx = s.dx
            phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
            s._frozenRho = s._face_rho(phix)
        x0 = np.concatenate([s.phi.ravel(), s.k.copy()])
        sol = newton_krylov(func, x0, iter=60, verbose=False,
                            f_tol=1e-5, f_rtol=1e-8)
        s.phi[:] = sol[:nphi].reshape(s.nx, s.ny, s.nz)
        s.k[:] = sol[nphi:]
        L = s.fp_loads() if fullpot else s.loads()
        print('  epoch %d: CL=%.4f k_mid=%.3f nsup=%d' % (
            epoch, L['CL'], s.k[len(s.k) // 2],
            int(s._sup_mask().sum())), flush=True)
    s._frozen = s._frozenA = s._frozenRho = None
    return s


if __name__ == '__main__':
    a = float(sys.argv[1]) if len(sys.argv) > 1 else 1.25
    run_m1(alpha_deg=a)
