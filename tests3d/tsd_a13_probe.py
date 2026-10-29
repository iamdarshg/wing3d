"""A13 probe: the mode-7 hybrid row, checked WITHOUT any TSD field solve.

Three cheap checks, all linear algebra on the 13x7 lattice:

  1. THE BIAS EXHIBIT. Put the exact 2D flat-plate doublet distribution on
     the lattice (mu = pi*alpha*F(x)) and read the downwash back two ways:
       D.mu  -- the mode-6 analytic kernel (must be -1*alpha: A12),
       E.mu  -- `w_plane`'s own estimator applied to the SAME sheet, i.e.
                the horseshoe field sampled at w_plane's two stencil
                heights (+0.0074c, -0.0246c). This is the estimator's
                systematic bias, measured as a chordwise curve.
     B = L.E - D is what mode 7 subtracts; B.mu is that bias in velocity
     units, so the hybrid row is bias-free by construction.

  2. THE CANCELLATION IS EXACT, NOT APPROXIMATE. Replace the field by its
     analytic equivalent (w_plane := wav = L.E.mu) and Newton-solve the
     mode-7 tangency system. The hybrid row collapses to D.mu - g, so the
     root must equal the mode-6 root to solver tolerance. This exercises
     the whole mode-7 assembly (E, L, B, the 3/4-chord map, the signs).

  3. THE COUPLING IS LIVE. d(row)/d(phi) by central difference on the four
     stencil nodes, against d(row)/d(mu), and the mode-6 control which must
     be identically zero.

Run: python tests3d/tsd_a13_probe.py
"""
import numpy as np
import sys

sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from wing3d.tsd import TSDSolver, flat_plate_jump          # noqa: E402
from scipy.optimize import newton_krylov                  # noqa: E402


def build(span, M=0.3, alpha_deg=1.25, ny=64, mode=7):
    s = TSDSolver(Minf=M, alpha_deg=alpha_deg, nx=80, ny=ny, nz=32,
                  span=span, linear=True)
    s.kutta_mode = mode
    s.kutta_on = True
    s.kutta_free_mu = True
    s.kutta_kernel = True
    s.kutta_hybrid = (mode == 7)
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_tangency = True
    s.kutta_te_kutta = True
    s.kutta_w = 1.0
    s.init_free_mu()
    s.dlm_kernel()
    return s


def gmat(s):
    nc, ng = s.n_ctrl, len(s.wing_jy)
    x_at = np.asarray(s._dlmInfo['xcp']) / s.chord
    sl = (np.asarray(s.profile['camber'](x_at)) if s.profile is not None
          else np.zeros(nc))
    return (sl - s.alpha)[:, None] * np.ones((1, ng))


def vstack_kutta(s, A, b):
    """Discrete Kutta row: mu(TE) - mu(TE-1) = 0 on the last ng rows."""
    nc, ng = s.n_ctrl, len(s.wing_jy)
    last = np.arange((nc - 1) * ng, nc * ng)
    A = A.copy()
    A[last, :] = 0.0
    A[last, last] = 1.0
    A[last, last - ng] = -1.0
    b = b.copy()
    b[last] = 0.0
    return A, b


def cl_circ(s):
    G = s.mu_b[-1]
    jy = np.asarray(s.wing_jy)
    return 2.0 * float(np.sum(G * s.dy[jy])) / (s.chord * s.span), G


def section1(span, tag):
    s = build(span)
    nc, ng = s.n_ctrl, len(s.wing_jy)
    a = s.alpha
    xc = np.asarray(s._dlmInfo['xcp'])
    mu = np.outer(flat_plate_jump(xc) * np.pi * a, np.ones(ng))
    D = s.dlm_kernel()
    E = s.dlm_field_equivalent()
    L = s.dlm_obs_matrix()
    B = s.dlm_hybrid_bias()
    wD = (D @ mu.ravel()).reshape(nc, ng)[:, ng // 2] / a
    wE = (L @ E @ mu.ravel()).reshape(nc, ng)[:, ng // 2] / a
    wB = (B @ mu.ravel()).reshape(nc, ng)[:, ng // 2] / a
    print('  %s: x_cp    ' % tag, ' '.join('%7.3f' % v for v in xc))
    print('  %s: D.mu/a  ' % tag, ' '.join('%7.2f' % v for v in wD))
    print('  %s: E.mu/a  ' % tag, ' '.join('%7.2f' % v for v in wE))
    print('  %s: B.mu/a  ' % tag, ' '.join('%7.2f' % v for v in wB))
    return s


def solve6(s):
    D = s.dlm_kernel()
    A, b = vstack_kutta(s, D, gmat(s).ravel())
    mu = np.linalg.solve(A, b).reshape(s.n_ctrl, len(s.wing_jy))
    s.mu_b = mu
    return mu


def solve7_ideal(s):
    """Mode-7 Newton solve with the field replaced by its analytic
    equivalent: w_plane := wav = L.E.mu, i.e. the exact doublet sheet."""
    nc, ng = s.n_ctrl, len(s.wing_jy)
    L = s.dlm_obs_matrix()
    E = s.dlm_field_equivalent()
    iu = np.asarray(s.wing_ix)
    jy = np.asarray(s.wing_jy)
    nn = len(iu)
    real_w = s.w_plane

    def fake_w_plane(wide=False):
        wnode = (E @ s.mu_b.ravel()).reshape(nn, ng)
        wav = np.zeros((s.nx, s.ny))
        wav[np.ix_(iu, jy)] = wnode
        return None, None, wav
    s.w_plane = fake_w_plane
    try:
        def f(y):
            s.mu_b = y.reshape(nc, ng)
            return s.hybrid_tangency_residual().ravel()
        y = newton_krylov(f, s.mu_b.ravel().copy(), iter=60, verbose=False,
                          f_tol=1e-12, f_rtol=1e-12)
        s.mu_b = y.reshape(nc, ng)
        res = float(np.abs(f(y)).max())
    finally:
        s.w_plane = real_w
    return s.mu_b, res


def main():
    print('== 1. bias exhibit: same exact flat-plate sheet through D and E ==')
    s40 = section1(40.0, 'span40')
    print('     (exact sheet: -1.00 everywhere; D is A12-validated;')
    print('      E is the mode-5 estimator geometry -> its own bias curve)')
    s16 = section1(1.6, 'AR1.6 ')

    print('\n== 2. mode-7 assembly reduces to mode 6 for an ideal field ==')
    for tag, span in (('span40', 40.0), ('AR1.6', 1.6)):
        sA = build(span)
        m6 = solve6(sA)
        sB = build(span)
        m7, r7 = solve7_ideal(sB)
        print('  %s: max|mu_hybrid - mu_kernel| = %.3e  (row |F|max %.1e)'
              % (tag, float(np.abs(m7 - m6).max()), r7))

    print('\n== 3. coupling diagnostics on the AR1.6 lattice ==')
    s = build(1.6)
    solve6(s)
    rng = np.random.default_rng(0)
    s.phi = 0.0 * s.phi
    rep0 = s.hybrid_coupling_report()
    print('  phi=0        : d(row)/d(phi)=%.3e d(row)/d(mu)=%.3e '
          'ratio=%.3f  mode6 dphi=%.3e'
          % (rep0['dphi_row'], rep0['dmu_row'], rep0['ratio_phi_mu'],
             rep0['dphi_row_mode6']))
    s.phi = rng.standard_normal(s.phi.shape) * 1e-3
    rep = s.hybrid_coupling_report()
    print('  phi=random   : d(row)/d(phi)=%.3e d(row)/d(mu)=%.3e '
          'ratio=%.3f  mode6 dphi=%.3e'
          % (rep['dphi_row'], rep['dmu_row'], rep['ratio_phi_mu'],
             rep['dphi_row_mode6']))
    print('  B: max|B|=%.4f  max row-sum|B|=%.4f  (D max=%.4f)'
          % (rep['B_max'], rep['B_l1_row_max'], rep['D_max']))

    print('\n== 4. mode-6 regression (must be unchanged by the mode-7 code) ==')
    a = np.radians(1.25)
    for tag, span in (('span40', 40.0), ('AR1.6', 1.6)):
        s2 = build(span)
        solve6(s2)
        CL, G = cl_circ(s2)
        print('  %s: CL_circ=%.5f  CL/alpha=%.3f rad^-1' %
              (tag, CL, CL / a))
    print('  (A12 recorded span40 0.852 of 2*pi*alpha, AR1.6 2.24/rad)')


if __name__ == '__main__':
    main()
