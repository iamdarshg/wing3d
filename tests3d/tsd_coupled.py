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


def build_case(Minf, alpha_deg, linear, fullpot=False, box=None,
               nx=80, ny=32, nz=32):
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
    s = TSDSolver(Minf=Minf, alpha_deg=alpha_deg, nx=nx, ny=ny, nz=nz,
                  omega=0.8, wake_gamma=g, linear=linear,
                  fullpot=fullpot, box=box)
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
    return s, res, mesh


def panel_aft_target(mesh, res, yc, xq=0.9):
    """Panel sectional dCp (lower-upper) at x=xq per span station yc.

    Physical Kutta target from the panel solution (no self-reference).
    """
    C = mesh.centroid
    cp = np.asarray(res['cp'])
    out = np.zeros(len(yc))
    for j, y in enumerate(yc):
        du, dl = None, None
        for sgn, tag in [(1, 'u'), (-1, 'l')]:
            m = (np.abs(C[:, 1] - y) < 0.35) & ((C[:, 2] * sgn) > 0.001)
            if not np.any(m):
                continue
            o = np.argsort(C[m, 0])
            xx, cc = C[m, 0][o], cp[m][o]
            v = float(np.interp(xq, xx, cc))
            if tag == 'u':
                du = v
            else:
                dl = v
        if du is not None and dl is not None:
            out[j] = dl - du
    return out


def run_coupled(Minf, alpha_deg=1.25, linear=True, kw=10.0, verbose=True,
                fullpot=False, epochs=10):
    s, res, mesh = build_case(Minf, alpha_deg, linear, fullpot)
    nphi = s.nx * s.ny * s.nz
    nk = len(s.wing_jy)
    s.k[:] = 1.0
    # panel aft-loading target (physical Kutta reference, no
    # self-reference): TSD must reproduce panel's sectional dCp at
    # x=0.9 while the pocket develops forward of it.
    tare = panel_aft_target(mesh, res, s.yc[s.wing_jy])
    if verbose:
        print('panel target mean=%.4f' % tare.mean(), flush=True)

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
        # Picard-outer (map+A/rho refresh) + NK-inner (frozen smooth epoch)
        for epoch in range(epochs):
            s._frozen = s._sup_mask()
            s._frozenA = s._face_coeff()
            if fullpot:
                dx = s.dx
                phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / \
                    dx[:, None, None]
                s._frozenRho = s._face_rho(phix)
            x0 = np.concatenate([s.phi.ravel(), s.k.copy()])
            sol = newton_krylov(func, x0, iter=60, verbose=False,
                                f_tol=1e-5, f_rtol=1e-8)
            s.phi[:] = sol[:nphi].reshape(s.nx, s.ny, s.nz)
            s.k[:] = sol[nphi:]
            L = s.fp_loads() if fullpot else s.loads()
            print('  epoch %d: CL=%.4f k_mid=%.3f nsup=%d' % (
                epoch, L['CL'], s.k[len(s.k) // 2],
                int(s._sup_mask().sum())), flush=True)
        s._frozen = None
        s._frozenA = None
        s._frozenRho = None
    L = s.fp_loads() if fullpot else s.loads()
    return s, L


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    linear = not (len(sys.argv) > 2 and sys.argv[2] == 'nl')
    fp = len(sys.argv) > 3 and sys.argv[3] == 'fp'
    a = 4.0 if linear else 1.25
    s, L = run_coupled(Minf, alpha_deg=a, linear=linear, fullpot=fp)
    print('coupled M%.2f %s%s: CL=%.4f CD=%.5f k_mid=%.3f' % (
        Minf, 'lin' if linear else 'NL', '/fp' if fp else '', L['CL'],
        L['CD'], s.k[len(s.k) // 2]))
