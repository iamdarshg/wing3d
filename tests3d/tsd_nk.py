"""TSD Newton-Krylov driver (robust alternative to SOR)."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                            wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.optimize import newton_krylov


def run_tsd(Minf, alpha_deg=4.0, nx=80, ny=32, nz=32, verbose=True):
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
                  omega=1.0, wake_gamma=g)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    s.wake_scale = 1.0
    # pin farfield box
    F = s.farfield
    s.phi[0, :, :] = F[0, :, :]
    s.phi[-1, :, :] = F[-1, :, :]
    s.phi[:, 0, :] = F[:, 0, :]
    s.phi[:, -1, :] = F[:, -1, :]
    s.phi[:, :, 0] = F[:, :, 0]
    s.phi[:, :, -1] = F[:, :, -1]
    N = s.nx * s.ny * s.nz

    def func(x):
        s.phi[:] = x.reshape(s.nx, s.ny, s.nz)
        R = s.residual()
        # freeze farfield (Dirichlet): zero residual there
        R[0, :, :] = 0
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        return R.ravel()

    x0 = s.phi.ravel().copy()
    sol = newton_krylov(func, x0, iter=100, verbose=verbose, f_tol=1e-6)
    s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
    L = s.loads()
    return L


if __name__ == '__main__':
    L = run_tsd(0.3)
    print('TSD M0.3 NK: CL=%.4f CD=%.5f (panel 0.367)' % (L['CL'], L['CD']))
