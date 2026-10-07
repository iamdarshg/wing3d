"""Picard(NK-inner) TSD: frozen map+A epochs solved exactly by Newton-Krylov."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.spatial import cKDTree
from scipy.optimize import newton_krylov


def build_case(Minf, alpha_deg=1.25):
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
                  omega=0.8, wake_gamma=g)
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
    return s


def func_factory(s):
    def func(x):
        s.phi[:] = x.reshape(s.nx, s.ny, s.nz)
        R = s.residual()
        R[0, :, :] = 0
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        return R.ravel()
    return func


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8
    s = build_case(Minf)
    s.wake_scale = 1.0
    func = func_factory(s)
    for epoch in range(12):
        # freeze map + full A (Picard linearization)
        s._frozen = s._sup_mask()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenA = s._coeff_A(phix)
        # Jacobi preconditioner (stretched-grid anisotropy ~1e6 cond)
        beta2 = max(1 - s.Minf ** 2, 0.05)
        ddx = np.append(s.dx, s.dx[-1])
        ddy = np.append(s.dy, s.dy[-1])
        ddz = np.append(s.dz, s.dz[-1])
        D = (2 * beta2 / ddx ** 2)[:, None, None] + \
            (2 / ddy ** 2)[None, :, None] + (2 / ddz ** 2)[None, None, :]
        from scipy.sparse.linalg import LinearOperator
        N = s.nx * s.ny * s.nz
        Minv = LinearOperator((N, N),
                              matvec=lambda v: (v / D.ravel()))
        x0 = s.phi.ravel().copy()
        sol = newton_krylov(func, x0, iter=200, verbose=False, f_tol=1e-4,
                            f_rtol=1e-10, inner_M=Minv)
        s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
        L = s.loads()
        nsup = int(s._sup_mask().sum())
        print('epoch %d: CL=%.4f CD=%.5f nsup=%d' % (epoch, L['CL'],
                                                    L['CD'], nsup),
              flush=True)
    s._frozen = None
    s._frozenA = None
    jmid = int(np.argmin(np.abs(s.yc)))
    print('midspan upper Cp(x):')
    for i in s.wing_ix[::4]:
        print('  x=%.3f Cp=%.3f' % (s.xc[i], L['cpu'][i, jmid]))
