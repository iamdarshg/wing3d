"""1D Kutta root-find: scale k on (wing jump + wake Gamma) to zero TE dCp."""
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
    base_jump = (mu_up[Tu.query(pts)[1]] -
                 mu_lo[Tl.query(pts)[1]]).reshape(XX.shape)
    base_G = s.wake_G.copy()
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    F = s.farfield
    s.phi[0, :, :] = F[0, :, :]
    s.phi[-1, :, :] = F[-1, :, :]
    s.phi[:, 0, :] = F[:, 0, :]
    s.phi[:, -1, :] = F[:, -1, :]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    return s, base_jump, base_G


def solve_epoch(s, nepoch=6):
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
    for _ in range(nepoch):
        if not s.linear:
            s._frozen = s._sup_mask()
            dx = s.dx
            phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
            s._frozenA = s._coeff_A(phix)
        x0 = s.phi.ravel().copy()
        sol = newton_krylov(func, x0, iter=60, verbose=False, f_tol=1e-4,
                            f_rtol=1e-10)
        s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
    s._frozen = None
    s._frozenA = None


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    linear = len(sys.argv) > 2 and sys.argv[2] == 'lin'
    s, base_jump, base_G = build_case(Minf, 1.25, linear)
    s.wake_scale = 1.0
    jmid = int(np.argmin(np.abs(s.yc)))
    ite = s.wing_ix[-1]

    def eval_k(k):
        s.wing_jump = k * base_jump
        s.wake_G = k * base_G
        s.phi[1:-1, 1:-1, 1:-1] = 0  # fresh start per eval
        F = s.farfield
        s.phi[0, :, :] = F[0, :, :]
        s.phi[-1, :, :] = F[-1, :, :]
        s.phi[:, 0, :] = F[:, 0, :]
        s.phi[:, -1, :] = F[:, -1, :]
        s.phi[:, :, 0] = F[:, :, 0]
        s.phi[:, :, -1] = F[:, :, -1]
        solve_epoch(s)
        L = s.loads()
        dcp = L['cpl'][ite, jmid] - L['cpu'][ite, jmid]
        return dcp, L['CL']

    # secant on f(k) = dcp_te(k); bracket [0.5, 1.5]
    k0, k1 = 0.5, 1.5
    f0, cl0 = eval_k(k0)
    print('k=%.2f dcp=%+.4f CL=%.4f' % (k0, f0, cl0), flush=True)
    f1, cl1 = eval_k(k1)
    print('k=%.2f dcp=%+.4f CL=%.4f' % (k1, f1, cl1), flush=True)
    for it in range(6):
        if abs(f1 - f0) < 1e-12:
            break
        k2 = float(np.clip(k1 - f1 * (k1 - k0) / (f1 - f0), 0.2, 4.0))
        f2, cl2 = eval_k(k2)
        print('k=%.3f dcp=%+.4f CL=%.4f' % (k2, f2, cl2), flush=True)
        k0, f0 = k1, f1
        k1, f1 = k2, f2
        if abs(f1) < 5e-3:
            break
