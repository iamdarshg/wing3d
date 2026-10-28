"""Measure the Kutta-row weight for the free-circulation TSD solve.

The coupled system is a saddle-point block::

    [ A   B ] [dphi]   [ R_phi ]
    [ C   0 ] [dGam] = [ R_kutt ]

A = dR_phi/dphi (the discrete TSD operator, magnitude ~ 1/dy^2+1/dz^2),
B = dR_phi/dGamma (the sheet source), C = dR_kutt/dphi (Cp is built
from phi; the Kutta rows have NO direct Gamma term). Eliminating dphi
gives the effective Kutta-row gain C A^-1 B ~ C*B/A, so the weight that
puts the Kutta rows on the same Jacobian footing as the phi rows is::

    w = A_scale^2 / (C_scale * B_scale)

Everything is measured by finite differences on the real g80x32x32 grid
(no arbitrary 10x). Usage: python tsd_kutta_weight.py [M] [geo]
"""
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

os.environ.setdefault('WING3D_IMPRINT', 'll')
os.environ.setdefault('WING3D_KUTTA', '1')


def main(M=0.3, geo='dfvlr_r4'):
    from tsd_campaign import build_profile, build_section
    from tsd_coupled import build_case, lifting_line_gamma
    from tsd_m_sup import run_sup  # noqa: F401  (import check)

    prof = build_profile(geo)
    sec = build_section(geo)
    s, _res, _mesh = build_case(M, 1.25, False, True, None, 80, 32, 32,
                                1.6, section=sec, profile=prof)
    s.Minf = M
    s.profile = prof
    s._frozen = s._frozenA = s._frozenRho = None
    s._frozenRhon = None
    s.kutta_mode = int(os.environ.get('WING3D_KUTTA', '1'))
    s.kutta_on = True
    s.kutta_bound = (s.kutta_mode == 1)
    s.bound_ramp = False
    s.wake_scale = 1.0
    jy = np.asarray(s.wing_jy)
    s.gamma_b = np.asarray(lifting_line_gamma(1.25, M, s.span,
                                             profile=prof, y=s.yc[jy]),
                           dtype=float)
    ng = len(jy)
    # a converged-ish phi helps the scales be representative: a few
    # line-relaxation sweeps on the seeded state (cheap, no NK)
    for _ in range(6):
        s.line_relax(omega=0.8)
    R0 = s.residual()
    K0 = s.kutta_residual(upstream=2)
    k0 = s.k0
    iu = np.asarray(s.wing_ix)
    i_te = iu[-3]

    # ---- A: phi-row diagonal (perturb a node, read its own row) ------
    d = 1e-6
    diags = []
    for (i, j) in [(i_te, jy[ng // 2]), (iu[2], jy[ng // 2]),
                   (i_te, jy[1]), (i_te, jy[-2])]:
        for k in (k0 + 1, k0 - 1):
            base = R0[i, j, k]
            s.phi[i, j, k] += d
            Rp = s.residual()
            s.phi[i, j, k] -= d
            diags.append((Rp[i, j, k] - base) / d)
    A_scale = float(np.median(np.abs(diags)))

    # ---- C: Kutta rows vs phi (same perturbation family) ------------
    # NOTE probe at the EXTRAPOLATION planes (k0+-1, k0+-2) and at
    # i_te+-1: the Kutta row is a central x-difference at i_te, so a
    # node at i_te itself does not appear in it (measured C=0 there).
    Cmax = 0.0
    for (i, j) in [(i_te + 1, jy[ng // 2]), (i_te - 1, jy[ng // 2]),
                   (i_te + 1, jy[1]), (i_te + 1, jy[-2])]:
        for k in (k0 + 1, k0 + 2, k0 - 1, k0 - 2):
            if not (0 < k < s.nz - 1):
                continue
            base = s.kutta_residual(upstream=2)
            s.phi[i, j, k] += d
            Kp = s.kutta_residual(upstream=2)
            s.phi[i, j, k] -= d
            Cmax = max(Cmax, float(np.abs((Kp - base) / d).max()))
    C_scale = Cmax

    # ---- B: phi rows vs Gamma (sheet source) ------------------------
    dG = 1e-5
    Bmax = 0.0
    Kg = 0.0
    for j in (ng // 2, 1, ng - 2):
        s.gamma_b[j] += dG
        Rp = s.residual()
        Kp = s.kutta_residual(upstream=2)
        s.gamma_b[j] -= dG
        Bmax = max(Bmax, float(np.abs((Rp - R0) / dG).max()))
        Kg = max(Kg, float(np.abs((Kp - K0) / dG).max()))
    B_scale = Bmax

    w = A_scale ** 2 / max(C_scale * B_scale, 1e-30)
    print('grid %dx%dx%d  dy_core=%.4f  dz_core=%.4f  dx_te=%.4f'
          % (s.nx, s.ny, s.nz, float(np.median(s.dy)), float(np.median(s.dz)),
             float(s.dx[i_te])))
    print('A_scale (dR_phi/dphi)  = %.4e' % A_scale)
    print('C_scale (dR_kutt/dphi) = %.4e' % C_scale)
    print('B_scale (dR_phi/dGam)  = %.4e' % B_scale)
    print('direct dRk/dGam (frozen phi, should be ~0) = %.2e' % Kg)
    print('effective Kutta gain C*B/A = %.4e' % (C_scale * B_scale / A_scale))
    # ---- effective coupled gain: perturb Gamma and CONVERGE phi ------
    # (two Newton-Krylov solves, so this is the true dR_kutt/dGamma
    # after the elliptic response, not the C*B/A upper bound)
    from scipy.optimize import newton_krylov
    nphi = (s.nx - 1) * s.ny * s.nz

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
        return s.kutta_residual(upstream=2)

    K0 = phi_solve()
    dG2 = 1e-4
    jm = ng // 2
    s.gamma_b[jm] += dG2
    K1 = phi_solve()
    s.gamma_b[jm] -= dG2
    geff = float((K1[jm] - K0[jm]) / dG2)
    print('converged coupled dCp(TE)/dG at mid = %.4e' % geff)
    print('linearised Kutta root: Gamma* = %.5f (LL mid %.5f)'
          % (s.gamma_b[jm] - K0[jm] / geff, s.gamma_b[jm]))
    w = A_scale / max(abs(geff), 1e-30)
    print('=> w = A_scale/|geff| = %.4g' % w)
    return w

    return w


if __name__ == '__main__':
    M = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    geo = sys.argv[2] if len(sys.argv) > 2 else 'dfvlr_r4'
    main(M, geo)
