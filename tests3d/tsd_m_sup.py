"""TSD full-potential at M>=1: sup inflow-pinned, outflow-free BC test.

Hypothesis from RCA: at M=1.0 the TSD stays on the subsonic branch
(CL ~ 0.115 incompressible-like) because ALL box faces are free-floating
early-branch DOFs. At M>=1 the flow is hyperbolic-in-x: the upstream
face must be pinned to uniform (phi=0), the downstream face outflow-free,
so the solution can develop the bow-shock/supersonic-pocket character
(near-symmetric, CL collapses, wave drag grows).

Runs M=1.0 and M=1.1 AR6, reports CL/net section + supersonic fraction.
Truth references (OF AR1.6 anchors, NOT AR-matched): M1.0 CL=0.0787,
M1.1 CL=0.0533. AR6 should sit ABOVE the AR1.6 anchor (tip relief).
"""
import numpy as np
import sys
import os
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_coupled import build_case
from scipy.optimize import newton_krylov


def run_sup(M, span=6.0, epochs=10, free_k=True, profile=None, section=None,
            base=None, callback=None, alpha_deg=None):
    kutta = int(os.environ.get('WING3D_KUTTA', '0') or 0)
    a_deg = float(os.environ.get('WING3D_ALPHA', '1.25')) \
        if alpha_deg is None else float(alpha_deg)
    if base is None:
        alpha_deg = a_deg
        s, res, mesh = build_case(M, alpha_deg, False, True, None, 80, 32,
                                  32, span, section=section,
                                  profile=profile)
    else:
        # reuse a prebuilt state: panel/imprint is Mach-independent, only
        # the operator coefficients depend on Minf (read live in residual)
        s = base
    s.Minf = M
    # clear stale frozen state: a previous Mach that raised leaves the
    # lagged mask/A/rho on the shared `base`, silently seeding this run
    # (perf-audit trap T5).
    s._frozen = s._frozenA = s._frozenRho = None
    s._frozenRhon = None
    if profile is not None:
        s.profile = profile
    # strip panel imprints that lock the subsonic circulation:
    # keep wing_jump (ramp shape) since free-k proved unstable, but
    # force new BCs
    # NOTE: s.farfield=None forces phi=0 on the whole outer box, which
    # forbids circulation (lift then comes only from transpiration and
    # cannot lose lift transonically). WING3D_FARFIELD=wake keeps the
    # wake-consistent lifting farfield computed by build_case. With a
    # FREE circulation (WING3D_KUTTA=1) boxing phi=0 fights the lift,
    # so the wake farfield is the default there.
    ff = os.environ.get('WING3D_FARFIELD', 'wake' if kutta else 'off')
    if ff == 'wake':
        s.phi[:] = 0.0
        F = s.farfield
        if F is not None:
            s.phi[0, :, :] = F[0, :, :]
            s.phi[-1, :, :] = F[-1, :, :]
            s.phi[:, 0, :] = F[:, 0, :]
            s.phi[:, -1, :] = F[:, -1, :]
            s.phi[:, :, 0] = F[:, :, 0]
            s.phi[:, :, -1] = F[:, :, -1]
    else:
        s.farfield = None
        s.phi[:] = 0.0
    s.k[:] = 1.0
    ng = len(s.wing_jy)
    if kutta:
        # FREE bound-vortex circulation: seed from lifting line (starting
        # at 0 or at the old panel value lands in the wrong basin).
        from tsd_coupled import lifting_line_gamma
        s.kutta_mode = kutta
        s.kutta_on = True
        s.kutta_bound = (kutta == 1)
        # the prescribed panel/LL bound sheet is replaced by the free one
        s.bound_ramp = False
        s.gamma_b = np.asarray(lifting_line_gamma(
            a_deg, M, s.span, chord=s.chord, profile=profile,
            y=s.yc[np.asarray(s.wing_jy)]), dtype=float)
        print('  kutta: Gamma0 (LL) mid=%.5f  2*sum(G*dy)/S=%.5f' % (
            s.gamma_b[len(s.gamma_b) // 2],
            2.0 * float(np.sum(s.gamma_b * s.dy[np.asarray(s.wing_jy)]))
            / (s.chord * s.span)), flush=True)
    else:
        s.kutta_on = False
        s.kutta_bound = False
    nx, ny, nz = s.nx, s.ny, s.nz
    nphi = (nx - 1) * ny * nz   # interior + outflow face (inflow excluded)
    nk = len(s.k)

    def unpack(y):
        phi = np.zeros((nx, ny, nz))
        phi[1:, :, :] = y[:nphi].reshape(nx - 1, ny, nz)
        # inflow face: uniform supersonic freestream (phi=0)
        phi[0, :, :] = 0.0
        return phi

    def func(y):
        s.phi[:] = unpack(y[:nphi])
        if kutta:
            s.gamma_b[:] = y[nphi:nphi + ng]
        elif free_k:
            s.k[:] = y[nphi:]
        R = s.residual()
        bc_rows = os.environ.get('WING3D_BC_ROWS', '') == '1'
        if bc_rows:
            # audit ablation: replace the rank-deficient zero rows with
            # real farfield statements (Dirichlet 0 on lateral/top/bottom,
            # linear extrapolation on the outflow plane).
            R[-1, :, :] = s.phi[-1, :, :] - s.phi[-2, :, :]
            R[:, 0, :] = s.phi[:, 0, :]
            R[:, -1, :] = s.phi[:, -1, :]
            R[:, :, 0] = s.phi[:, :, 0]
            R[:, :, -1] = s.phi[:, :, -1]
        else:
            # out-flow face free (zero residual row), lateral/top/bottom free
            R[-1, :, :] = 0
            R[:, 0, :] = 0
            R[:, -1, :] = 0
            R[:, :, 0] = 0
            R[:, :, -1] = 0
        R[0, :, :] = 0   # inflow pinned by unpack, no equation needed
        if kutta:
            # circulation rows: TE dCp(j) = 0, weighted to the phi-row
            # Jacobian scale (measured, see tsd_kutta_weight.py)
            Rk = s.kutta_residual(upstream=2) * s.kutta_w
        elif free_k:
            Rk = s.kutta_residual(upstream=2) * 10.0
        else:
            Rk = None
        # only return the free DOF rows (inflow plane excluded)
        if kutta or free_k:
            return np.concatenate([R[1:].ravel(), Rk])
        return R[1:].ravel()

    s.kutta_hist = []
    for ep in range(epochs):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenRho = s._face_rho(phix)
        s._frozenRhon = s._nodal_rho()
        if kutta:
            y0 = np.concatenate([s.phi[1:].ravel(), s.gamma_b])
        elif free_k:
            y0 = np.concatenate([s.phi[1:].ravel(), s.k.copy()])
        else:
            y0 = s.phi[1:].ravel()
        y = newton_krylov(func, y0, iter=80, verbose=False,
                          f_tol=1e-4, f_rtol=1e-7)
        if kutta:
            s.phi[:] = unpack(y[:nphi])
            s.gamma_b[:] = y[nphi:nphi + ng]
            resn = float(np.abs(func(y)).max())
        elif free_k:
            s.phi[:] = unpack(y[:nphi])
            s.k[:] = y[nphi:]
            resn = float('nan')
        else:
            s.phi[:] = unpack(y)
            resn = float('nan')
        L = s.fp_loads()
        supn = int(s._sup_mask().sum())
        if kutta:
            gmid = s.gamma_b[len(s.gamma_b) // 2]
            kres = float(np.abs(s.kutta_residual(upstream=2)).max())
            cl_circ = 2.0 * float(np.sum(
                s.gamma_b * s.dy[np.asarray(s.wing_jy)])) / (
                    s.chord * s.span)
            s.kutta_hist.append({'ep': ep, 'CL': float(L['CL']),
                                 'gamma_mid': float(gmid),
                                 'CL_circ': cl_circ,
                                 'kutta_max': kres, 'res_max': resn,
                                 'nsup': supn})
            print('  ep %d: CL=%.4f (2$G dy/S=%.4f) G_mid=%.4f '
                  'Kmax=%.2e |F|max=%.1e nsup=%d' % (
                      ep, L['CL'], cl_circ, gmid, kres, resn, supn),
                  flush=True)
        else:
            print('  ep %d: CL=%.4f k_mid=%.3f nsup=%d' % (
                ep, L['CL'], s.k[len(s.k) // 2], supn), flush=True)
        if callback is not None:
            callback(ep, float(L['CL']), supn)
        # viscous-inviscid coupling (WING3D_BL=1): after each inviscid
        # epoch, re-march the integral BL on the current surface Cp and
        # feed the relaxed displacement back through the wing
        # transpiration; the next epoch re-solves the inviscid problem.
        if getattr(s, 'bl_on', False):
            s.bl_update(verbose=True)
    s._frozen = s._frozenA = s._frozenRho = None
    s._frozenRhon = None
    L = s.fp_loads()
    if getattr(s, 'bl_on', False):
        L['bl'] = {'hist': s.bl_hist, 'sep': s.bl_sep,
                   'dstar_u': s.bl_dstar_u, 'dstar_l': s.bl_dstar_l}
    if kutta:
        L['kutta'] = {'hist': s.kutta_hist, 'Gamma': s.gamma_b.copy(),
                      'w': s.kutta_w,
                      'y': s.yc[np.asarray(s.wing_jy)].copy()}
    jy = np.asarray(s.wing_jy)
    j0 = jy[np.argmin(np.abs(s.yc[jy]))]
    ix = np.asarray(s.wing_ix)
    xc = s.xc[ix] / s.chord
    cl = float(np.trapezoid((L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc))
    print('M=%.2f: full-wing CL=%.4f  centerline cl=%.4f  Cfmax=%.2f' % (
        M, L['CL'], cl, L['cpu'][ix, j0].max()), flush=True)
    return L, xc, ix, j0


if __name__ == '__main__':
    for M in [1.0, 1.1]:
        print('=== M=%.2f (AR6, inflow-pinned/outflow-free) ===' % M,
              flush=True)
        run_sup(M)