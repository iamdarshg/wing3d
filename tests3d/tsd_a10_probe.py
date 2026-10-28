"""A10 probe: does a TANGENCY row carry information about gamma_b?

A9 established that the old closure (dCp at 87% chord) and the Morino-style
wake condition are both structurally dead. This measures, on the real
g80x32x32 AR1.6 grid at M0.3, the two candidate rows of the new closure:

  TANGENCY   wav(x_cp, y_j) - (slope_c - alpha) = 0
             wav = 0.5*(wu + wl) = the field's normal velocity at the mean
             plane (tsd.TSDSolver.w_plane). Reported for both transpiration
             variants: mode 1 (thickness only -> the lifting normal velocity
             is NOT imposed anywhere, so the row is live) and mode 3
             (alpha+camber in the transpiration -> the row is imposed by the
             fixed source and can carry no information).

  MORINO     gamma_b(j) - [phi_up(TE) - phi_lo(TE)] = 0
             measured as te_jump(gamma) vs gamma. If dJ/dgamma == 1 the row
             is an identity in this discretization (the cut face flux IS the
             definition of the jump) and cannot select the amplitude.

For each convergent phi (sheet amplitude frozen) it reports CL, CL_circ,
the span profile of wav, and the finite-difference gains dwav/dgamma and
dJ/dgamma, which are the numbers the row weight is built from.

Usage: python tsd_a10_probe.py [M] [geo]
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))
os.environ.setdefault('WING3D_IMPRINT', 'll')


def main(M=0.3, geo='dfvlr_r4'):
    from scipy.optimize import newton_krylov
    from tsd_campaign import build_profile, build_section
    from tsd_coupled import build_case, lifting_line_gamma

    sec = build_section(geo)
    s, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6,
                           section=sec, profile=None)
    s.Minf = M
    s.bound_ramp = False
    s.wake_scale = 1.0
    s.wake_sign = 1.0
    jy = np.asarray(s.wing_jy)
    ng = len(jy)
    G0 = np.asarray(lifting_line_gamma(1.25, M, s.span, profile=None,
                                       y=s.yc[jy]), dtype=float)
    s.gamma_b = G0.copy()
    # NOTE: the sheet source only exists when the mode is on; set mode 1
    # (thickness-only transpiration + bound sheet) BEFORE measuring the
    # Jacobian scales, else dR_phi/dGamma measures 0.
    s.kutta_mode = 1
    s.kutta_on = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    # Jacobian scales of the phi block (same measurement as
    # tsd_kutta_weight.py: perturb a node, read its own row)
    s._frozen = s._sup_mask()
    s._frozenA = s._face_coeff()
    phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / s.dx[:, None, None]
    s._frozenRho = s._face_rho(phix)
    s._frozenRhon = s._nodal_rho()
    R0 = s.residual()
    k0 = s.k0
    iu = np.asarray(s.wing_ix)
    i_te = iu[-3]
    d = 1e-6
    diags = []
    for (i, j) in [(i_te, jy[ng // 2]), (iu[2], jy[ng // 2]),
                   (i_te, jy[1]), (i_te, jy[-2])]:
        for k in (k0 + 1, k0 - 1):
            base = R0[i, j, k]
            s.phi[i, j, k] += d
            diags.append((s.residual()[i, j, k] - base) / d)
            s.phi[i, j, k] -= d
    A_scale = float(np.median(np.abs(diags)))
    # B: phi rows vs gamma (sheet source), frozen phi
    dG = 1e-5
    Bmax = 0.0
    for p in (ng // 2, 1, ng - 2):
        s.gamma_b[p] += dG
        Bmax = max(Bmax, float(np.abs((s.residual() - R0) / dG).max()))
        s.gamma_b[p] -= dG
    s._frozen = s._frozenA = s._frozenRho = None
    s._frozenRhon = None
    print('A_scale (dR_phi/dphi) = %.4e   B_scale (dR_phi/dG) = %.4e'
          % (A_scale, Bmax), flush=True)

    jmid = ng // 2
    jm = int(jy[jmid])          # grid column of the mid-span wing cell
    # control station: wing cell nearest 3/4 chord
    xw = s.xc[iu] / s.chord
    i_cp = int(iu[np.argmin(np.abs(xw - 0.75))])
    print('wing_ix %d cells, x/c = %s' % (len(iu), np.round(xw, 3)),
          flush=True)
    print('control station i_cp = %d  x/c = %.3f' % (i_cp, s.xc[i_cp]),
          flush=True)

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
        return float(np.abs(f(y)).max())

    def measure(mode, scale):
        s.kutta_mode = mode
        s.kutta_on = True
        s.kutta_bound = mode in (1, 3)
        s.kutta_thick_only = (mode == 1)
        s.gamma_b = scale * G0
        t0 = time.time()
        res = phi_solve()
        wu, wl, wav = s.w_plane()
        L = s.fp_loads()
        J = s.te_jump()
        # geometric mean-plane slope at the stations
        xc = s.xc[iu] / s.chord
        sc_ = np.zeros((len(iu), 1))
        if s.profile is not None:
            sc_ = s.profile['camber'](xc)[:, None]
        geom = sc_ - s.alpha
        cl_circ = 2.0 * float(np.sum(s.gamma_b * s.dy[jy])) / (
            s.chord * s.span)
        row = {'mode': mode, 'scale': float(scale),
               'CL': float(L['CL']), 'CL_circ': cl_circ,
               'G_mid': float(s.gamma_b[jmid]),
               'res': res, 'secs': time.time() - t0,
               'wav_mid': float(wav[i_cp, jm]),
               'wav_span': [float(wav[i_cp, j]) for j in jy],
               'wav_chord_mid': [float(wav[i, jm]) for i in iu],
               'geom_cp': float(geom[np.argmin(np.abs(xc - 0.75)), 0]),
               'J_mid': float(J[jm]),
               'J_span': [float(J[jy[p]]) for p in range(ng)],
               'Gspan': [float(s.gamma_b[p]) for p in range(ng)]}
        print('  mode %d s=%+.2f: CL=%+.4f CLcirc=%+.4f G_mid=%+.4f  '
              'wav_cp_mid=%+.5f  J_mid=%+.5f  (res %.1e, %.0fs)' % (
                  mode, scale, row['CL'], cl_circ, row['G_mid'],
                  row['wav_mid'], row['J_mid'], res, row['secs']),
              flush=True)
        print('     wav(x/c) mid: %s' % np.array2string(
            np.array(row['wav_chord_mid']), precision=4), flush=True)
        return row

    rows = []
    rows.append(measure(1, 0.5))
    rows.append(measure(1, 1.0))
    rows.append(measure(1, 1.5))
    rows.append(measure(3, 0.0))
    rows.append(measure(3, 1.0))

    # ---- frozen-phi gains of the two candidate rows -------------------
    s.kutta_mode = 1
    s.kutta_thick_only = True
    s.kutta_bound = True
    s.gamma_b = G0.copy()
    phi_solve()
    wu, wl, wav = s.w_plane()
    J = s.te_jump()
    tot = 0.0
    for p in (jmid, 1, ng - 2):
        j = int(jy[p])
        s.gamma_b[p] += dG
        wu2, wl2, wav2 = s.w_plane()
        J2 = s.te_jump()
        s.gamma_b[p] -= dG
        tot = max(tot, float(abs((wav2[i_cp, j] - wav[i_cp, j]) / dG)))
        print('  frozen-phi d wav/dG (cp,col=%d) = %+.4e   d J/dG = %+.4e'
              % (j, float((wav2[i_cp, j] - wav[i_cp, j]) / dG),
                 float((J2[j] - J[j]) / dG)), flush=True)
    print('  max frozen-phi |d wav/dG| at cp = %.4e' % tot, flush=True)

    out = {'M': M, 'geo': geo, 'A_scale': A_scale, 'B_scale': Bmax,
           'i_cp': int(i_cp), 'x_cp': float(s.xc[i_cp]),
           'xw': xw.tolist(), 'y': s.yc[jy].tolist(), 'rows': rows}
    p = os.path.join(REPO, 'tests3d', 'a10_probe_M%.3f_%s.json' % (M, geo))
    json.dump(out, open(p, 'w'), indent=1)
    print('saved', p)


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.3,
         sys.argv[2] if len(sys.argv) > 2 else 'dfvlr_r4')
