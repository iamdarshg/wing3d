"""A10 probe 3: spanwise shape and gain of the tangency row.

The tangency row (mode 4) is  wav(x_cp, y_j) - (dyc/dx - alpha) = 0  with
wav the SYMMETRIC first-order estimator of the plane normal velocity
(tsd.TSDSolver.w_plane, which cancels the +-J/2 jump imprint that makes
the wide extrapolation read 24.4*gamma).

Two things decide the coupled answer and both are measured here, at fixed
sheet amplitude (phi converged):
  * the spanwise shape of wav at fixed gamma. If it is flat, the row is
    span-diagonal and the selected gamma must come out span-flat (A9's
    failure mode); if it falls at the tips the row couples the span and
    elliptic loading can emerge.
  * the per-station gain d wav/d gamma (converged, i.e. phi re-solved),
    which is the number the row weight is built from.

It also prints the *linear* prediction of the fixed point
gamma_j* = (w_geom_j - wav0_j) / gain_j next to the lifting-line value.

Usage: python tsd_a10_probe3.py [M]
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))
os.environ.setdefault('WING3D_IMPRINT', 'll')


def main(M=0.3, geo='dfvlr_r4'):
    from scipy.optimize import newton_krylov
    from tsd_campaign import build_profile, build_section
    from tsd_coupled import build_case, lifting_line_gamma

    sec = build_section(geo) if geo != 'sym' else None
    prof = build_profile(geo) if geo != 'sym' else None
    s, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6,
                           section=sec, profile=prof)
    if prof is not None:
        s.profile = prof
    s.Minf = M
    s.bound_ramp = False
    s.wake_scale = 1.0
    s.kutta_mode = 4
    s.kutta_on = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_tangency = True
    s.kutta_morino = True
    s.tang_xc = 0.75
    jy = np.asarray(s.wing_jy)
    ng = len(jy)
    G0 = np.asarray(lifting_line_gamma(1.25, M, s.span, profile=None,
                                       y=s.yc[jy]), dtype=float)
    iu = np.asarray(s.wing_ix)
    xw = s.xc[iu] / s.chord
    i_cp = int(iu[int(np.argmin(np.abs(xw - 0.75)))])
    # geometric mean-plane slope requirement at the control station
    if s.profile is not None:
        w_geom = float(s.profile['camber'](s.xc[i_cp] / s.chord) - s.alpha)
    else:
        w_geom = float(-s.alpha)

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

    out = {'M': M, 'x_cp': float(s.xc[i_cp]), 'y': s.yc[jy].tolist(),
           'w_geom': float(w_geom), 'G0': G0.tolist(), 'scales': {}}
    for sc in (0.0, 1.0, 1.5):
        s.gamma_b = sc * G0
        phi_solve()
        wu, wl, wav = s.w_plane()
        s.tang_mode = 'pt'
        r_pt = np.asarray(s.tangency_residual())
        s.tang_mode = 'gal'
        r_gal = np.asarray(s.tangency_residual())
        wav_j = np.array([wav[i_cp, j] for j in jy])
        L = s.fp_loads()
        out['scales']['%.2f' % sc] = {
            'wav': wav_j.tolist(), 'CL': float(L['CL']),
            'r_pt': r_pt.tolist(), 'r_gal': r_gal.tolist(),
            'tau_row': float(np.abs(r_pt).max())}
        print('s=%.2f CL=%.4f  wav_j = %s' % (
            sc, L['CL'], np.array2string(wav_j, precision=4)), flush=True)
        print('     row pt  = %s' % np.array2string(
            r_pt, precision=4), flush=True)
        print('     row gal = %s' % np.array2string(
            r_gal, precision=4), flush=True)
    w0 = np.array(out['scales']['0.00']['wav'])
    w1 = np.array(out['scales']['1.00']['wav'])
    w15 = np.array(out['scales']['1.50']['wav'])
    gain = (w15 - w0) / (1.5 * G0)
    gain2 = (w1 - w0) / G0
    print('per-station converged gain d wav/dG (s0->1.0): %s' % (
        np.array2string(gain2, precision=3)), flush=True)
    print('per-station converged gain d wav/dG (s0->1.5): %s' % (
        np.array2string(gain, precision=3)), flush=True)
    g_pt = (np.array(out['scales']['1.00']['r_pt']) -
            np.array(out['scales']['0.00']['r_pt'])) / G0
    g_gal = (np.array(out['scales']['1.00']['r_gal']) -
             np.array(out['scales']['0.00']['r_gal'])) / G0
    r_gal0 = np.array(out['scales']['0.00']['r_gal'])
    print('ROW pt  converged gain dR/dG = %s' % (
        np.array2string(g_pt, precision=3)), flush=True)
    print('ROW gal converged gain dR/dG = %s' % (
        np.array2string(g_gal, precision=3)), flush=True)
    print('ROW gal root  G* = R0/(-gain) = %s' % (np.array2string(
        -r_gal0 / np.where(np.abs(g_gal) > 1e-12, g_gal, np.nan),
        precision=4), ), flush=True)
    pred = (w_geom - w0) / np.where(np.abs(gain2) > 1e-9, gain2, np.nan)
    print('linear tangency root G* = %s' % (
        np.array2string(pred, precision=4)), flush=True)
    print('lifting line G0        = %s' % (
        np.array2string(G0, precision=4)), flush=True)
    print('w_geom (needed wav)    = %.4f' % w_geom, flush=True)
    out['gain'] = gain.tolist()
    out['gain2'] = gain2.tolist()
    out['G_pred'] = pred.tolist()
    json.dump(out, open(os.path.join(
        REPO, 'tests3d', 'a10_probe3_M%.3f.json' % M), 'w'), indent=1)
    print('A_scale-consistent weight w = 7.24e3 / median|gain| = %.4g'
          % (7.24e3 / max(float(np.median(np.abs(gain2))), 1e-12)),
          flush=True)
    print('A_scale-consistent weight for the gal row = 7.24e3 / median|g_gal|'
          ' = %.4g' % (7.24e3 / max(float(np.median(np.abs(g_gal))), 1e-12)),
          flush=True)


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.3,
         sys.argv[2] if len(sys.argv) > 2 else 'dfvlr_r4')
