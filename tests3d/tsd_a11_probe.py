"""A11 probe: MEASURED row scaling for the mode-5 chordwise tangency system.

The mode-5 closure is one tangency row per (chordwise panel, span cell):
    row[i,j] = wav(x_i, y_j) - (dyc/dx(x_i) - alpha)
Rows read phi only, so (as in A9/A10) the frozen-phi derivative
d row / d mu is exactly 0 and the honest measure is the CONVERGED
sensitivity: re-solve phi at two doublet matrices and difference.

  python tsd_a11_probe.py sym      uncambered (profile=None), M0.3, AR1.6
  python tsd_a11_probe.py cam      dfvlr_r4 profile, M0.3, AR1.6

Reports, per geometry:
  * A_scale = median |dR_phi/dphi| (the phi-row diagonal scale the
    closure rows have to be weighted against)
  * D[i,j]  = converged d wav / d mu per control station and span cell
    (secant between mu = 0.6 and mu = 1.4 x the LL seed F(x)*Gamma_LL)
  * the resulting scalar weight w = A_scale/median|D| and the spread
Nothing is tuned here: the weight is derived from the measured D.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import newton_krylov

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_campaign import build_profile, build_section   # noqa: E402
from tsd_coupled import build_case, lifting_line_gamma  # noqa: E402
from wing3d.tsd import flat_plate_jump                  # noqa: E402
from tsd_m_sup import run_sup                           # noqa: E402


def setup(M=0.3, span=1.6, profile=None, section='dfvlr_r4'):
    sec = build_section(section)
    s, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                           span, section=sec, profile=profile)
    s.Minf = M
    s._frozen = s._frozenA = s._frozenRho = s._frozenRhon = None
    s.profile = profile
    s.k[:] = 1.0
    # mode-5 closure state (same flags run_sup sets)
    s.kutta_mode = 5
    s.kutta_on = True
    s.kutta_free_mu = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_tangency = True
    s.bound_ramp = False
    s.init_free_mu()
    return s


def phi_solve(s, mu, ep=2):
    """Converged phi at a FIXED doublet matrix mu (mode-5 rows inactive)."""
    s.mu_b[:] = mu
    nx, ny, nz = s.nx, s.ny, s.nz
    nphi = (nx - 1) * ny * nz

    def unpack(y):
        phi = np.zeros((nx, ny, nz))
        phi[1:, :, :] = y.reshape(nx - 1, ny, nz)
        return phi

    def func(y):
        s.phi[:] = unpack(y)
        R = s.residual()
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        R[0, :, :] = 0
        return R[1:].ravel()

    for _e in range(ep):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenRho = s._face_rho(phix)
        s._frozenRhon = s._nodal_rho()
        y = newton_krylov(func, s.phi[1:].ravel(), iter=80, verbose=False,
                          f_tol=1e-4, f_rtol=1e-7)
        s.phi[:] = unpack(y)
    return s


def a_scale(s, nprobe=8):
    """median |dR_phi/dphi| from one-sided phi differences at wing nodes."""
    k0 = s.k0
    iu = np.asarray(s.wing_ix)
    jy = np.asarray(s.wing_jy)
    base = s.residual().copy()
    vals = []
    eps = 1e-6
    for i in np.unique(np.linspace(0, len(iu) - 1, nprobe).astype(int)):
        for kk in (k0, k0 - 1):
            for j in (jy[0], jy[len(jy) // 2]):
                idx = (int(iu[i]), int(j), int(kk))
                s.phi[idx] += eps
                R = s.residual()
                s.phi[idx] -= eps
                d = (R[idx] - base[idx]) / eps
                if abs(d) > 0:
                    vals.append(abs(d))
    return float(np.median(vals)), vals


def probe(kind='sym'):
    prof = None if kind == 'sym' else build_profile('dfvlr_r4')
    t0 = time.time()
    s = setup(0.3, 1.6, profile=prof)
    a_deg = 1.25
    g_ll = np.asarray(lifting_line_gamma(
        a_deg, 0.3, s.span, chord=s.chord, profile=prof,
        y=s.yc[np.asarray(s.wing_jy)]), dtype=float)
    F = flat_plate_jump(s.x_ctrl)
    nc, ng = s.n_ctrl, len(s.wing_jy)
    mu_seed = F[:, None] * g_ll[None, :]
    print('kind=%s nc=%d ng=%d x_ctrl=%s' % (
        kind, nc, ng, np.round(s.x_ctrl, 3).tolist()), flush=True)
    print('  Gamma_LL mid=%.5f  mu_TE mid=%.5f' % (
        g_ll[len(g_ll) // 2], mu_seed[-1, len(g_ll) // 2]), flush=True)

    out = {}
    A_sc, vals = a_scale(s)
    out['A_scale'] = A_sc
    out['A_scale_vals'] = vals
    print('  A_scale (median |dR/dphi|) = %.3e  (n=%d)' % (A_sc, len(vals)),
          flush=True)

    iu = np.asarray(s.wing_ix)[s.i_ctrl][:, None]
    jj = np.asarray(s.wing_jy)[None, :]
    phi0 = s.phi.copy()
    wavs = {}
    for amp, tag in ((0.6, 'lo'), (1.4, 'hi')):
        s.phi[:] = phi0
        t = time.time()
        phi_solve(s, amp * mu_seed, ep=2)
        _wu, _wl, wav = s.w_plane()
        wavs[tag] = wav[iu, jj]
        print('  solved phi at amp=%.1f (%.0fs)  wav mid=%.5f' % (
            amp, time.time() - t, wavs[tag][nc // 2, ng // 2]), flush=True)

    dmu = 0.8 * mu_seed
    D = (wavs['hi'] - wavs['lo']) / dmu
    out['D'] = D.tolist()
    out['wav_lo'] = wavs['lo'].tolist()
    out['x_ctrl'] = s.x_ctrl.tolist()
    out['y'] = s.yc[np.asarray(s.wing_jy)].tolist()
    out['Gamma_LL'] = g_ll.tolist()
    out['mu_seed'] = mu_seed.tolist()
    med = float(np.median(np.abs(D)))
    out['D_median'] = med
    out['w_scalar'] = A_sc / med if med > 0 else None
    out['secs'] = time.time() - t0
    print('  D (converged d wav/d mu), rows = x_ctrl, cols = span cells:')
    for i in range(nc):
        print('   x=%.3f  ' % s.x_ctrl[i] +
              ' '.join('%7.2f' % D[i, j] for j in range(ng)), flush=True)
    print('  median|D|=%.3f  min=%.3f max=%.3f  -> w = A/|D| = %.3e' % (
        med, np.abs(D).min(), np.abs(D).max(), A_sc / med), flush=True)
    path = os.path.join(REPO, 'tests3d', 'a11_probe_%s.json' % kind)
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print('  -> %s (%.0fs)' % (path, out['secs']), flush=True)
    return out


if __name__ == '__main__':
    probe(sys.argv[1] if len(sys.argv) > 1 else 'sym')
