"""A10 probe 2: what does dphi/dz near the wing plane actually contain?

tsd_a10_probe.py found wav(x_cp) = 24.4 * gamma at M0.3 (mid-span). That
is 1/dz_core: the value is dominated by the potential JUMP carried on the
cut, not by the induced downwash (which lifting line says is O(0.01) for
AR1.6, CL 0.059 => alpha_i = CL/(pi AR) = 0.012).

This isolates it: dump the z-profile of phi at the control station on the
wing, at a station downstream (wake only) and at a station outboard of
the tip (no sheet at all), and compare four estimators of w at the plane:

  w_extrap  = (phi[k0+2]-phi[k0]) / (z[k0+2]-z[k0])  (tsd.w_plane)
  w_1s      = (phi[k0+1]-phi[k0]) / (z[k0+1]-z[k0])  (first-order, above)
  w_through = (phi[k0]-phi[k0-1] - jump) / dz        (mass through the cut)
  w_far     = (phi[k0+4]-phi[k0+3]) / dz             (4 cells up, no sheet)

and prints them for gamma = 0, 0.5, 1.0 Gamma_LL so the gamma-response of
each estimator is visible. That response, not the raw value, is what a
tangency row would have to use.

Usage: python tsd_a10_probe2.py [M]
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))
os.environ.setdefault('WING3D_IMPRINT', 'll')


def main(M=0.3):
    from scipy.optimize import newton_krylov
    from tsd_campaign import build_section
    from tsd_coupled import build_case, lifting_line_gamma
    from wing3d.tsd import flat_plate_jump

    sec = build_section('dfvlr_r4')
    s, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6,
                           section=sec, profile=None)
    s.Minf = M
    s.bound_ramp = False
    s.wake_scale = 1.0
    s.kutta_mode = 1
    s.kutta_on = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_tangency = False
    jy = np.asarray(s.wing_jy)
    ng = len(jy)
    jm = int(jy[ng // 2])
    G0 = np.asarray(lifting_line_gamma(1.25, M, s.span, profile=None,
                                       y=s.yc[jy]), dtype=float)
    s.gamma_b = G0.copy()
    k0 = s.k0
    iu = np.asarray(s.wing_ix)
    xw = s.xc[iu] / s.chord
    i_cp = int(iu[np.argmin(np.abs(xw - 0.75))])
    # a station just downstream of the TE (wake only) and one outboard of
    # the tip (no sheet at all)
    i_wake = min(int(s.wake_ix[2]), s.nx - 2) if len(s.wake_ix) else i_cp
    j_out = int(np.where(np.abs(s.yc) > 0.5 * s.span + 0.35)[0][0])
    print('k0=%d  z[k0-2..k0+2]=%s  dz[k0-2..k0+1]=%s' % (
        k0, np.array2string(s.z[max(k0 - 2, 0):k0 + 3], precision=4),
        np.array2string(s.dz[max(k0 - 2, 0):k0 + 2], precision=4)), flush=True)
    print('i_cp=%d (x/c %.3f)  i_wake=%d  outboard col=%d (y=%.2f)'
          % (i_cp, s.xc[i_cp], i_wake, j_out, s.yc[j_out]), flush=True)

    def phi_solve():
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / s.dx[:, None, None]
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

        y = newton_krylov(f, s.phi[1:].ravel(), iter=80, verbose=False,
                          f_tol=1e-7, f_rtol=1e-8)
        f(y)
        s._frozen = s._frozenA = s._frozenRho = None
        s._frozenRhon = None

    rows = []
    for sc in (0.0, 0.5, 1.0):
        s.gamma_b = sc * G0
        phi_solve()
        jump = flat_plate_jump(s.xc[i_cp] / s.chord) * s.gamma_b[ng // 2]
        dphi = s.phi[i_cp, jm, k0] - s.phi[i_cp, jm, k0 - 1]
        w_extrap = (s.phi[i_cp, jm, k0 + 2] - s.phi[i_cp, jm, k0]) / (
            s.z[k0 + 2] - s.z[k0])
        w_1s = (s.phi[i_cp, jm, k0 + 1] - s.phi[i_cp, jm, k0]) / (
            s.z[k0 + 1] - s.z[k0])
        w_thr = (dphi - jump) / s.dz[k0 - 1]
        w_far = (s.phi[i_cp, jm, k0 + 4] - s.phi[i_cp, jm, k0 + 3]) / (
            s.z[k0 + 4] - s.z[k0 + 3])
        # symmetric first-order pair straddling the cut: the +-J/2 imprint
        # of the potential jump enters with OPPOSITE sign on the two sides,
        # so the average cancels it to first order -- this is the estimator
        # that can see the induced downwash (the others are swamped by
        # jump/dz, see a10_probe_report)
        w_sym = 0.5 * ((s.phi[i_cp, jm, k0 + 1] - s.phi[i_cp, jm, k0]) / (
            s.z[k0 + 1] - s.z[k0]) + (s.phi[i_cp, jm, k0 - 1]
                                      - s.phi[i_cp, jm, k0 - 2]) / (
            s.z[k0 - 1] - s.z[k0 - 2]))
        prof_w = [float(s.phi[i_cp, jm, k]) for k in range(k0 - 3, k0 + 4)]
        prof_k = [float(s.phi[i_wake, jm, k]) for k in range(k0 - 3, k0 + 4)]
        prof_o = [float(s.phi[i_cp, j_out, k]) for k in range(k0 - 3, k0 + 4)]
        L = s.fp_loads()
        row = {'scale': sc, 'G': float(s.gamma_b[ng // 2]),
               'CL': float(L['CL']), 'jump': float(jump),
               'dphi': float(dphi), 'w_extrap': float(w_extrap),
               'w_1s': float(w_1s), 'w_through': float(w_thr),
               'w_sym': float(w_sym), 'w_far': float(w_far),
               'prof_wing': prof_w,
               'prof_wake': prof_k, 'prof_out': prof_o}
        rows.append(row)
        print('s=%.2f G=%+.4f jump=%+.4f dphi=%+.4f | w_extrap=%+.4f '
              'w_1s=%+.4f w_sym=%+.5f w_through=%+.5f w_far=%+.5f CL=%.4f'
              % (sc, row['G'], row['jump'], row['dphi'], w_extrap, w_1s,
                 w_sym, w_thr, w_far, row['CL']), flush=True)
        print('   z-prof wing : %s' % np.array2string(
            np.array(prof_w), precision=4), flush=True)
        print('   z-prof wake : %s' % np.array2string(
            np.array(prof_k), precision=4), flush=True)
        print('   z-prof outbd : %s' % np.array2string(
            np.array(prof_o), precision=4), flush=True)
    out = {'M': M, 'k0': int(k0), 'rows': rows,
           'z': s.z[max(k0 - 4, 0):k0 + 5].tolist()}
    p = os.path.join(REPO, 'tests3d', 'a10_probe2_M%.3f.json' % M)
    json.dump(out, open(p, 'w'), indent=1)
    print('saved', p)


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.3)
