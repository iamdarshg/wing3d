"""Gamma-scan of the Kutta row: is there a physical root and a basin?

For a fixed scale s, set gamma_b = s * Gamma_LL, converge phi with the
sheet amplitude FROZEN, and report CL, the TE loading dCp(TE) and the
circulation lift 2*sum(g dy)/S. This maps the closure function
f(Gamma) = dCp(TE) whose zero the coupled solve looks for, and shows
whether the lifting-line amplitude is a root at all.

Usage: python tsd_kutta_scan.py [M] [mode] [scales...]
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))


def main(M=0.3, mode=3, scales=(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)):
    os.environ['WING3D_KUTTA'] = str(mode)
    os.environ.setdefault('WING3D_IMPRINT', 'll')
    from scipy.optimize import newton_krylov
    from tsd_campaign import build_section
    from tsd_coupled import build_case, lifting_line_gamma
    sec = build_section('dfvlr_r4')
    s, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6,
                           section=sec, profile=None)
    s.Minf = M
    s.kutta_mode = mode
    s.kutta_on = True
    s.kutta_bound = mode in (1, 3)
    s.kutta_thick_only = (mode == 1)
    s.bound_ramp = False
    s.wake_scale = 1.0
    s._frozen = s._frozenA = s._frozenRho = None
    s._frozenRhon = None
    jy = np.asarray(s.wing_jy)
    G0 = np.asarray(lifting_line_gamma(1.25, M, s.span, profile=None,
                                       y=s.yc[jy]), dtype=float)
    s.gamma_b = G0.copy()

    def phi_solve():
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenRho = s._face_rho(phix)
        s._frozenRhon = s._nodal_rho()

        def f(y):
            s.phi[1:, :, :] = y.reshape(s.nx - 1, s.ny, s.nz)
            s.phi[0, :, :] = 0.0
            R = s.residual()
            R[-1, :, :] = 0
            R[:, 0, :] = 0
            R[:, -1, :] = 0
            R[:, :, 0] = 0
            R[:, :, -1] = 0
            R[0, :, :] = 0
            return R[1:].ravel()

        y = newton_krylov(f, s.phi[1:].ravel(), iter=60, verbose=False,
                          f_tol=1e-7, f_rtol=1e-8)
        f(y)
        s._frozen = s._frozenA = s._frozenRho = None
        s._frozenRhon = None

    rows = []
    for sc in scales:
        s.gamma_b = sc * G0
        phi_solve()
        K = s.kutta_residual(upstream=2)
        L = s.fp_loads()
        cl_circ = 2.0 * float(np.sum(s.gamma_b * s.dy[jy])) / (s.chord * s.span)
        row = {'scale': float(sc), 'CL': float(L['CL']), 'CL_circ': cl_circ,
               'dCpTE_mid': float(K[len(K) // 2]), 'dCpTE_max': float(K.max()),
               'dCpTE_absmean': float(np.abs(K).mean())}
        rows.append(row)
        print('  s=%+.2f  G_mid=%+.4f  CL=%+.4f  CLcirc=%+.4f  '
              'dCpTE mid=%+.4f  |mean|=%.4f' % (
                  sc, s.gamma_b[len(s.gamma_b) // 2], row['CL'],
                  cl_circ, row['dCpTE_mid'], row['dCpTE_absmean']),
              flush=True)
    p = os.path.join(REPO, 'tests3d', 'a9_gscan_M%.3f_m%d.json' % (M, mode))
    json.dump(rows, open(p, 'w'), indent=1)
    print('saved', p)


if __name__ == '__main__':
    M = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    mode = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    sc = tuple(float(v) for v in sys.argv[3:]) or \
        (-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
    main(M, mode, sc)
