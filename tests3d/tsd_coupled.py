"""Coupled TSD: Newton-Krylov over [phi, k] with Kutta rows.

k = per-span-cell bound-vortex amplitudes (free); Kutta rows
(TE dCp = 0) close the system. Validates linear (k -> 1) then
runs nonlinear transonic.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.spatial import cKDTree
from scipy.optimize import newton_krylov


def build_case(Minf, alpha_deg, linear):
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(alpha_deg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    meta = mesh.meta
    nl, ns = meta['n_loop'], meta['n_span']
    strips, _ = wing_strips_from_structured(nl, ns)
    info = [(up, lo, meta['te_points'][j], meta['te_points'][j + 1])
            for (up, lo, j) in strips]
    wakes_long = build_wake(mesh, info, length=9.0, n_panels=6)
    g = {'y': np.array([mesh.centroid[j * 48][1] for j in range(24)]),
         'gamma': np.array(res['muw'])}
    s = TSDSolver(Minf=Minf, alpha_deg=alpha_deg, nx=80, ny=32, nz=32,
                  omega=0.8, wake_gamma=g, linear=linear)
    C = mesh.centroid
    mu = np.asarray(res['mu'])
    up_idx = C[:, 2] > 0.001
    lo_idx = C[:, 2] < -0.001
    Tu, Tl = cKDTree(C[up_idx][:, :2]), cKDTree(C[lo_idx][:, :2])
    mu_up, mu_lo = mu[up_idx], mu[lo_idx]
    xc = s.xc[s.wing_ix]
    yc = s.yc[s.wing_jy]
    XX, YY = np.meshgrid(xc, yc, indexing='ij')
    pts = np.stack([XX.ravel(), YY.ravel()], axis=-1)
    s.wing_jump = (mu_up[Tu.query(pts)[1]] -
                   mu_lo[Tl.query(pts)[1]]).reshape(XX.shape)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    F = s.farfield
    s.phi[0, :, :] = F[0, :, :]
    s.phi[-1, :, :] = F[-1, :, :]
    s.phi[:, 0, :] = F[:, 0, :]
    s.phi[:, -1, :] = F[:, -1, :]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    s.wake_scale = 1.0
    return s


def run_coupled(Minf, alpha_deg=1.25, linear=True, kw=10.0, verbose=True):
    s = build_case(Minf, alpha_deg, linear)
    nphi = s.nx * s.ny * s.nz
    nk = len(s.wing_jy)
    s.k[:] = 1.0
    # phase 0 (tare): linear prescribed solution's TE pattern.
    # In linear mode this IS the answer (k must stay 1); in nonlinear
    # mode the rows drive k to hold this pattern (Kutta up to the
    # shared panel discretization defect).
    s_lin = build_case(Minf, alpha_deg, True)
    s_lin.wing_jump = s.wing_jump
    s_lin.wake_scale = 1.0

    def func_lin(x):
        s_lin.phi[:] = x.reshape(s_lin.nx, s_lin.ny, s_lin.nz)
        R = s_lin.residual()
        R[0, :, :] = 0
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        return R.ravel()

    sol = newton_krylov(func_lin, s_lin.phi.ravel().copy(), iter=60,
                        verbose=False, f_tol=1e-4, f_rtol=1e-10)
    s_lin.phi[:] = sol.reshape(s_lin.nx, s_lin.ny, s_lin.nz)
    tare = s_lin.kutta_residual(upstream=2)
    if verbose:
        print('tare mean=%.4f' % tare.mean(), flush=True)

    def func(x):
        s.phi[:] = x[:nphi].reshape(s.nx, s.ny, s.nz)
        s.k[:] = x[nphi:]
        R = s.residual()
        R[0, :, :] = 0
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        Rk = (s.kutta_residual(upstream=2) - tare) * kw
        return np.concatenate([R.ravel(), Rk])

    x0 = np.concatenate([s.phi.ravel(), s.k.copy()])
    if linear:
        sol = newton_krylov(func, x0, iter=120, verbose=verbose, f_tol=1e-6,
                            f_rtol=1e-6)
        s.phi[:] = sol[:nphi].reshape(s.nx, s.ny, s.nz)
        s.k[:] = sol[nphi:]
    else:
        # Picard-outer (map+A refresh) + NK-inner (frozen smooth epoch)
        for epoch in range(10):
            s._frozen = s._sup_mask()
            dx = s.dx
            phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
            s._frozenA = s._coeff_A(phix)
            x0 = np.concatenate([s.phi.ravel(), s.k.copy()])
            sol = newton_krylov(func, x0, iter=60, verbose=False,
                                f_tol=1e-5, f_rtol=1e-8)
            s.phi[:] = sol[:nphi].reshape(s.nx, s.ny, s.nz)
            s.k[:] = sol[nphi:]
            L = s.loads()
            print('  epoch %d: CL=%.4f k_mid=%.3f nsup=%d' % (
                epoch, L['CL'], s.k[len(s.k) // 2],
                int(s._sup_mask().sum())), flush=True)
        s._frozen = None
        s._frozenA = None
    L = s.loads()
    return s, L


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    linear = not (len(sys.argv) > 2 and sys.argv[2] == 'nl')
    a = 4.0 if linear else 1.25
    s, L = run_coupled(Minf, alpha_deg=a, linear=linear)
    print('coupled M%.2f %s: CL=%.4f CD=%.5f k_mid=%.3f' % (
        Minf, 'lin' if linear else 'NL', L['CL'], L['CD'],
        s.k[len(s.k) // 2]))
