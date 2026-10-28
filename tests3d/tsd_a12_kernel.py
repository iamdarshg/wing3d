"""A12: validate the mode-6 analytic doublet-lattice kernel WITHOUT a TSD solve.

Section 1 -- the decisive validation: put the EXACT 2D flat-plate doublet
distribution on the lattice (mu(x) = Gamma*F(x), Gamma = pi*alpha) and read
the downwash back through D. On a large span the answer must be
w = -alpha everywhere (that is the definition of the flat-plate solution,
and it is what the mode-5 field estimator failed to reproduce: A11 measured
0.45*alpha at the LE to 1.7*alpha aft).

Section 2 -- matrix structure (self, neighbours, row sum, conditioning).

Section 3 -- solve D.mu = -(alpha) with the discrete Kutta row:
  * span 40 (2D limit): CL must -> 2*pi*alpha
  * AR1.6: compare with Prandtl lifting line (lifting_line_gamma)

Run: python tests3d/tsd_a12_kernel.py
"""
import numpy as np
import sys

sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from wing3d.tsd import TSDSolver, flat_plate_jump
from tsd_coupled import lifting_line_gamma


def build(span, M=0.3, alpha_deg=1.25, ny=64):
    s = TSDSolver(Minf=M, alpha_deg=alpha_deg, nx=80, ny=ny, nz=32,
                  span=span, linear=True)
    s.kutta_mode = 6
    s.kutta_on = True
    s.kutta_free_mu = True
    s.kutta_kernel = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_tangency = True
    s.kutta_te_kutta = True
    s.init_free_mu()
    return s


def solve_mu(s):
    """Linear solve of the mode-6 tangency system (Kutta row included)."""
    D = s.dlm_kernel()
    nc, ng = s.n_ctrl, len(s.wing_jy)
    info = s._dlmInfo
    x_at = np.asarray(info['xcp']) / s.chord
    sl = (np.asarray(s.profile['camber'](x_at)) if s.profile is not None
          else np.zeros(nc))
    g = (sl - s.alpha)[:, None] * np.ones((1, ng))
    A = D.copy()
    b = g.ravel().copy()
    last = np.arange((nc - 1) * ng, nc * ng)
    A[last, :] = 0.0
    A[last, last] = 1.0
    A[last, last - ng] = -1.0
    b[last] = 0.0
    mu = np.linalg.solve(A, b)
    s.mu_b = mu.reshape(nc, ng)
    return s.mu_b, np.linalg.cond(A)


def cl_circ(s):
    # mode 6: mu(x=0)=0 by construction (the LE cell ramps 0 -> mu_0), so
    # the section circulation is mu(TE) itself, not mu(TE)-mu(first station)
    G = s.mu_b[-1]
    jy = np.asarray(s.wing_jy)
    return 2.0 * float(np.sum(G * s.dy[jy])) / (s.chord * s.span), G


if __name__ == '__main__':
    a = np.radians(1.25)

    print('== 1. exact flat-plate distribution through D (span 40) ==')
    s0 = build(40.0)
    D = s0.dlm_kernel()
    nc, ng = s0.n_ctrl, len(s0.wing_jy)
    ex = s0._dlmInfo['ex']
    xc = s0._dlmInfo['xcp']
    mu = np.outer(flat_plate_jump(xc) * np.pi * a, np.ones(ng))
    w = (D @ mu.ravel()).reshape(nc, ng)
    jm = ng // 2
    print('   x_cp   ', ' '.join('%7.3f' % v for v in xc))
    print('   w/alpha', ' '.join('%7.2f' % v for v in w[:, jm] / a))
    print('   (exact doublet sheet: -1 everywhere; A11 field estimator:'
          ' -0.45 LE .. -1.7 aft)')

    print('\n== 2. structure (AR1.6, 13x7 lattice) ==')
    s1 = build(1.6)
    D1 = s1.dlm_kernel()
    nc, ng = s1.n_ctrl, len(s1.wing_jy)
    jm = (nc // 2) * ng + ng // 2
    print('   D chord neighbours= %s' % np.array2string(
        np.array([D1[jm + d * ng, jm] for d in (-2, -1, 1, 2)]), precision=4))
    print('   D span  neighbours= %s' % np.array2string(
        np.array([D1[jm, jm + d] for d in (-1, 1)]), precision=4))
    print('   cond(D_incl_Kutta)= %.3e'
          % np.linalg.cond(np.vstack([D1[:-ng], np.eye(nc * ng)[-ng:]
                                      - np.eye(nc * ng, k=-ng)[-ng:]])))

    print('\n== 3a. 2D limit (span 40): CL -> 2*pi*alpha = %.5f =='
          % (2 * np.pi * a))
    s2 = build(40.0)
    mu2, c2 = solve_mu(s2)
    CL2, G2 = cl_circ(s2)
    print('   CL_circ = %.5f  ratio = %.4f  cond = %.2e'
          % (CL2, CL2 / (2 * np.pi * a), c2))

    print('\n== 3b. AR1.6 (the ladder case) ==')
    s3 = build(1.6)
    mu3, c3 = solve_mu(s3)
    CL3, G3 = cl_circ(s3)
    yg = np.asarray(s3.yc[np.asarray(s3.wing_jy)])
    gam = lifting_line_gamma(1.25, 0.3, 1.6, chord=1.0, y=yg)
    CLll = 2.0 * float(np.sum(gam * s3.dy[np.asarray(s3.wing_jy)])) / \
        (s3.chord * s3.span)
    g_sec = np.diff(mu3, axis=0) / np.diff(np.asarray(s3.x_ctrl))[:, None]
    jmid = len(yg) // 2
    print('   mu rows x_cp  ', ' '.join('%7.3f' % v
                                        for v in s3._dlmInfo['xcp']))
    print('   gamma (per c) ', ' '.join('%7.4f' % v for v in g_sec[:, jmid]))
    print('   G/Gmid        ', ' '.join('%7.3f' % v for v in G3 / G3[jmid]))
    print('   LL G/Gmid     ', ' '.join('%7.3f' % v for v in gam / gam[jmid]))
    print('   kernel CL_circ = %.5f | LL CL = %.5f (%.1f%%) | cond %.2e'
          % (CL3, CLll, 100 * (CL3 / CLll - 1), c3))
