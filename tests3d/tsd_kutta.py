"""TSD with Kutta-updated wake: Gamma(y) from TE jump each epoch."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.optimize import newton_krylov


def build_case(Minf, alpha_deg=1.25, linear=True):
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
    # SEED wake with panel-linear Gamma (breaks symmetry); Kutta update
    # lets it grow to the nonlinear value
    g0 = {'y': np.array([mesh.centroid[j * 48][1] for j in range(24)]),
          'gamma': np.array(res['muw'])}
    s = TSDSolver(Minf=Minf, alpha_deg=alpha_deg, nx=80, ny=32, nz=32,
                  omega=0.8, wake_gamma=g0, linear=linear,
                  bound_ramp=False)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    F = s.farfield
    # start farfield at zero wake (ramp will grow); keep box shape via F*scale
    s.phi[0, :, :] = 0
    s.phi[-1, :, :] = 0
    s.phi[:, 0, :] = 0
    s.phi[:, -1, :] = 0
    s.phi[:, :, 0] = 0
    s.phi[:, :, -1] = 0
    return s, res


def te_jump(s):
    """Potential jump across wing plane at TE column, per y-cell."""
    k0 = s.k0
    # TE x-cell: last wing_ix
    i = s.wing_ix[-1]
    d = s.phi[i, :, k0 + 1] - s.phi[i, :, k0 - 2]
    return d


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
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    linear = len(sys.argv) > 2 and sys.argv[2] == 'lin'
    s, res = build_case(Minf, linear=linear)
    s.wake_scale = 1.0
    func = func_factory(s)
    relax = 0.5
    Kp = 2.0
    jmid = int(np.argmin(np.abs(s.yc)))
    i90 = s.wing_ix[np.argmin(np.abs(s.xc[s.wing_ix] - 0.9))]
    for epoch in range(15):
        s._frozen = None if linear else s._sup_mask()
        if not linear:
            dx = s.dx
            phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
            s._frozenA = s._coeff_A(phix)
        x0 = s.phi.ravel().copy()
        sol = newton_krylov(func, x0, iter=40, verbose=False, f_tol=1e-5,
                            f_rtol=1e-4)
        s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
        # pressure-Kutta sensor: aft loading -> grow/shrink Gamma
        # (measured at 90% chord; TE cell polluted by cut singularity)
        L = s.loads()
        dcp_90 = L['cpl'][i90, jmid] - L['cpu'][i90, jmid]
        fac = float(np.clip(1 + relax * Kp * dcp_90, 0.5, 2.0))
        s.wake_G = s.wake_G * fac
        print('epoch %d: CL=%.4f dcp90=%+.4f fac=%.3f maxG=%.4f'
              % (epoch, L['CL'], dcp_90, fac, np.abs(s.wake_G).max()),
              flush=True)
    s._frozen = None
    s._frozenA = None
