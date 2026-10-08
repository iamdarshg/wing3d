"""TSFOIL-style Kutta loop: NK-inner epochs + RECIRC (extrapolated TE jump)."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
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
                  omega=0.8, wake_gamma=g, linear=linear,
                  bound_ramp=False)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    F = s.farfield
    s.phi[0, :, :] = F[0, :, :]
    s.phi[-1, :, :] = F[-1, :, :]
    s.phi[:, 0, :] = F[:, 0, :]
    s.phi[:, -1, :] = F[:, -1, :]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    s.wake_scale = 1.0
    # farfield circulation target = panel TE Gamma per span cell
    s._circff = s.wake_G[s.wing_jy] if len(s.wing_jy) else None
    return s


def run_recirc(Minf, alpha_deg=1.25, linear=True, epochs=12):
    s = build_case(Minf, alpha_deg, linear)

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

    for ep in range(epochs):
        if not linear:
            s._frozen = s._sup_mask()
            dx = s.dx
            phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
            s._frozenA = s._coeff_A(phix)
        x0 = s.phi.ravel().copy()
        # Jacobi preconditioner (cut-dipole stiffness stalls plain NK)
        beta2 = max(1 - s.Minf ** 2, 0.05)
        ddx = np.append(s.dx, s.dx[-1])
        ddy = np.append(s.dy, s.dy[-1])
        ddz = np.append(s.dz, s.dz[-1])
        D = (2 * beta2 / ddx ** 2)[:, None, None] + \
            (2 / ddy ** 2)[None, :, None] + (2 / ddz ** 2)[None, None, :]
        from scipy.sparse.linalg import LinearOperator
        N = s.nx * s.ny * s.nz
        Minv = LinearOperator((N, N), matvec=lambda v: v / D.ravel())
        sol = newton_krylov(func, x0, iter=200, verbose=False, f_tol=1e-6,
                            f_rtol=1e-10, inner_M=Minv)
        s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
        dg = s.recirc_update(W=0.25)
        L = s.loads()
        print('ep %d: CL=%.4f dG=%.4f maxG=%.4f' % (
            ep, L['CL'], dg, np.abs(s.wake_G).max()), flush=True)
    s._frozen = None
    s._frozenA = None
    return s


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    linear = not (len(sys.argv) > 2 and sys.argv[2] == 'nl')
    a = 4.0 if linear else 1.25
    run_recirc(Minf, alpha_deg=a, linear=linear)
