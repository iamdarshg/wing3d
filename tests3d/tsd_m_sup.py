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
from wing3d.tsd import flat_plate_jump
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
        if not kutta:
            s.phi[:] = 0.0
        # kutta mode KEEPS the build_case interior (the lifting-line
        # imprint state, consistent with Gamma0): the task's own warning
        # is that starting from phi=0 lands in the wrong basin, and it
        # does -- measured, Gamma walked from +0.041 to -0.009.
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
    free_mu = (kutta in (5, 6, 7))
    if kutta:
        # FREE bound-vortex circulation: seed from lifting line (starting
        # at 0 or at the old panel value lands in the wrong basin).
        from tsd_coupled import lifting_line_gamma
        s.kutta_mode = kutta
        s.kutta_on = True
        s.kutta_free_mu = free_mu
        s.kutta_kernel = (kutta in (6, 7))
        s.kutta_hybrid = (kutta == 7)
        s.kutta_bound = kutta in (1, 3, 4, 5, 6, 7)
        s.kutta_thick_only = kutta in (1, 4, 5, 6, 7)
        s.kutta_tangency = kutta in (4, 5, 6, 7)
        s.kutta_morino = (kutta == 4) and (
            os.environ.get('WING3D_MORINO', '1') == '1')
        s.kutta_te_kutta = (kutta in (5, 6, 7)) and (
            os.environ.get('WING3D_KUTTA_TE', '1') == '1')
        s.kutta_wte = float(os.environ.get(
            'WING3D_KUTTA_WTE', '1.0' if kutta in (6, 7) else '6.1'))
        s.tang_xc = float(os.environ.get('WING3D_TANG_XC', '0.75'))
        s.tang_mode = os.environ.get('WING3D_TANG_MODE', 'gal')
        # the prescribed panel/LL bound sheet is replaced by the free one
        s.bound_ramp = False
        g_ll = np.asarray(lifting_line_gamma(
            a_deg, M, s.span, chord=s.chord, profile=profile,
            y=s.yc[np.asarray(s.wing_jy)]), dtype=float)
        s.gamma_b = g_ll
        if free_mu:
            # seed mu(x,y) with the cumulative flat-plate shape at the
            # lifting-line amplitude: mu = F(x) * Gamma_LL(y)  (the same
            # starting point mode 4 uses, spread over the chord)
            s.init_free_mu()
            s.mu_b = (flat_plate_jump(s.x_ctrl)[:, None]
                      * g_ll[None, :]).astype(float)
            print('  kutta5: mu seed F(x)*LL, nc=%d ng=%d, mu_TE mid=%.5f'
                  % (s.n_ctrl, len(s.wing_jy), s.mu_b[-1, len(g_ll) // 2]),
                  flush=True)
        print('  kutta: Gamma0 (LL) mid=%.5f  2*sum(G*dy)/S=%.5f' % (
            g_ll[len(g_ll) // 2],
            2.0 * float(np.sum(g_ll * s.dy[np.asarray(s.wing_jy)]))
            / (s.chord * s.span)), flush=True)
    else:
        s.kutta_on = False
        s.kutta_bound = False
        s.kutta_free_mu = False
    nx, ny, nz = s.nx, s.ny, s.nz
    nphi = (nx - 1) * ny * nz   # interior + outflow face (inflow excluded)
    nmu = int(s.mu_b.size) if free_mu else 0
    if free_mu:
        ng_tot = nmu
    else:
        ng_tot = ng
    nk = len(s.k)

    def unpack(y):
        phi = np.zeros((nx, ny, nz))
        phi[1:, :, :] = y[:nphi].reshape(nx - 1, ny, nz)
        # inflow face: uniform supersonic freestream (phi=0)
        phi[0, :, :] = 0.0
        return phi

    def tang_rows(st):
        """The mu-block rows of the active closure (mode 5/6/7)."""
        if getattr(st, 'kutta_hybrid', False):
            return st.hybrid_tangency_residual()
        if getattr(st, 'kutta_kernel', False):
            return st.dlm_tangency_residual()
        return st.tangency_panel_residual()

    def func(y):
        s.phi[:] = unpack(y[:nphi])
        if free_mu:
            s.mu_b[:] = y[nphi:nphi + nmu].reshape(s.n_ctrl, ng)
        elif kutta:
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
            # circulation rows: mode 1-3 TE dCp(j) = 0, mode 4 the
            # TANGENCY row wav(x_cp) - (dyc/dx - alpha) per span cell,
            # mode 5 one TANGENCY ROW PER (chordwise panel, span cell)
            # (all weighted to the phi-row Jacobian scale, measured --
            # see tests3d/tsd_kutta_weight.py, tsd_a10_probe.py and
            # tsd_a11_probe.py)
            if free_mu:
                if getattr(s, 'kutta_hybrid', False):
                    # mode 7: hybrid row w_field(phi) - B.mu - g. Linear in
                    # mu, phi-DEPENDENT (the Jacobian is no longer
                    # block-triangular), and the estimator's chordwise bias
                    # is removed by B.
                    Rk = s.hybrid_tangency_residual().ravel() * s.kutta_w
                elif getattr(s, 'kutta_kernel', False):
                    # mode 6: the rows are the ANALYTIC doublet-lattice
                    # tangency equation, and the discrete Kutta row is
                    # already row nc-1 of the same array (in mu units).
                    Rk = s.dlm_tangency_residual().ravel() * s.kutta_w
                else:
                    Rk = s.tangency_panel_residual().ravel() * s.kutta_w
                    if s.kutta_te_kutta:
                        # last chordwise station = the TE: replace its
                        # tangency row by the Kutta condition gamma(TE)=0,
                        # i.e. mu(TE) - mu(last-1) = 0, in mu units
                        # weighted by the MEASURED local self-response
                        _i0 = s.kutta_w * s.kutta_wte
                        Rk[-ng:] = (s.mu_b[-1, :] - s.mu_b[-2, :]) * _i0
            elif getattr(s, 'kutta_tangency', False):
                Rk = s.tangency_residual() * s.kutta_w
            else:
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
            if free_mu:
                y0 = np.concatenate([s.phi[1:].ravel(), s.mu_b.ravel()])
            else:
                y0 = np.concatenate([s.phi[1:].ravel(), s.gamma_b])
        elif free_k:
            y0 = np.concatenate([s.phi[1:].ravel(), s.k.copy()])
        else:
            y0 = s.phi[1:].ravel()
        y = newton_krylov(func, y0, iter=80, verbose=False,
                          f_tol=1e-4, f_rtol=1e-7)
        if free_mu:
            s.phi[:] = unpack(y[:nphi])
            s.mu_b[:] = y[nphi:nphi + nmu].reshape(s.n_ctrl, ng)
            resn = float(np.abs(func(y)).max())
        elif kutta:
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
            if free_mu:
                Gsec = s.mu_circulation()
                gmid = float(Gsec[len(Gsec) // 2])
                if getattr(s, 'kutta_hybrid', False) or \
                        getattr(s, 'kutta_kernel', False):
                    _tr = tang_rows(s)
                    kres = float(np.abs(_tr[:-1]).max())
                    kte = float(np.abs(_tr[-1]).max())
                else:
                    _tr = s.tangency_panel_residual()
                    kres = float(np.abs(_tr[:-1] if s.kutta_te_kutta
                                        else _tr).max())
                    kte = float(np.abs(s.mu_b[-1, :] - s.mu_b[-2, :]).max())
                _J = s.te_jump()
                jm = int(np.asarray(s.wing_jy)[len(s.wing_jy) // 2])
                morino = float(_J[jm] - s.mu_b[-1, len(s.gamma_b) // 2])
                extra = ' tang=%.2e KuttaTE=%.2e Morino d=%.2e muLE=%.4f' % (
                    kres, kte, morino, s.mu_b[0, len(s.gamma_b) // 2])
                cl_circ = 2.0 * float(np.sum(
                    Gsec * s.dy[np.asarray(s.wing_jy)])) / (
                        s.chord * s.span)
            else:
                gmid = s.gamma_b[len(s.gamma_b) // 2]
                if getattr(s, 'kutta_tangency', False):
                    kres = float(np.abs(s.tangency_residual()).max())
                    _J = s.te_jump()
                    _jy = np.asarray(s.wing_jy)
                    jm = int(_jy[len(_jy) // 2])
                    morino = float(_J[jm] - s.gamma_b[len(s.gamma_b) // 2])
                    extra = ' tang=%.2e Morino d=%.2e' % (kres, morino)
                else:
                    kres = float(np.abs(s.kutta_residual(upstream=2)).max())
                    morino = float('nan')
                    extra = ''
                cl_circ = 2.0 * float(np.sum(
                    s.gamma_b * s.dy[np.asarray(s.wing_jy)])) / (
                        s.chord * s.span)
            s.kutta_hist.append({'ep': ep, 'CL': float(L['CL']),
                                 'gamma_mid': float(gmid),
                                 'CL_circ': cl_circ,
                                 'kutta_max': kres, 'res_max': resn,
                                 'morino_mid': morino,
                                 'kutta_te': (float(kte) if free_mu else None),
                                 'nsup': supn})
            print('  ep %d: CL=%.4f (2$G dy/S=%.4f) G_mid=%.4f '
                  'Kmax=%.2e |F|max=%.1e nsup=%d%s' % (
                      ep, L['CL'], cl_circ, gmid, kres, resn, supn, extra),
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
        if free_mu:
            s.gamma_b = s.mu_circulation()
        L['kutta'] = {'hist': s.kutta_hist, 'Gamma': s.gamma_b.copy(),
                      'w': s.kutta_w,
                      'y': s.yc[np.asarray(s.wing_jy)].copy()}
        if free_mu:
            L['kutta']['mu'] = s.mu_b.copy()
            L['kutta']['x_ctrl'] = s.x_ctrl.copy()
            L['kutta']['Gamma_mu'] = s.mu_circulation().copy()
            L['kutta']['loading'] = s.mu_loading().copy()
            if getattr(s, 'kutta_hybrid', False):
                # A13: the row decomposition, so the bias cancellation and
                # the phi loop can be judged from the record instead of a
                # rerun. w_field = the field estimator, w_eq = its analytic
                # equivalent L.E.mu, bias = B.mu (what mode 7 subtracts).
                _wu, _wl, _wav = s.w_plane()
                _iu = np.asarray(s.wing_ix)
                _jy = np.asarray(s.wing_jy)
                _nc, _ng = s.n_ctrl, len(_jy)
                L['kutta']['w_field'] = np.asarray(
                    _wav[np.ix_(_iu, _jy)], dtype=float).copy()
                L['kutta']['w_eq'] = (s.dlm_field_equivalent()
                                      @ s.mu_b.ravel()).reshape(_nc, _ng)
                L['kutta']['bias'] = (s.dlm_hybrid_bias()
                                      @ s.mu_b.ravel()).reshape(_nc, _ng)
                L['kutta']['coupling'] = s.hybrid_coupling_report()
            L['kutta']['tang_rows'] = (tang_rows(s).copy()
                                       if (getattr(s, 'kutta_kernel', False)
                                           or getattr(s, 'kutta_hybrid',
                                                      False))
                                       else s.tangency_panel_residual().copy())
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