"""A14: the "1.36x" measured DIRECTLY, without a coupled solve.

A13 inferred the discrete field's normalisation from a converged mode-7
residual (`resid_phys`), which mixes three things: the thickness response,
the sheet response and iteration noise (|F|max was still 1e-2 at ep4).
This probe separates them: it PRESCRIBES an exact doublet sheet
(mu = pi*alpha*F(x), the flat-plate distribution mode 6's kernel was
validated on), solves the phi problem ONCE with that sheet frozen, and
reads

    ratio = < w_field / (E.mu) >      (mid-span, |E.mu| > 0.25 alpha)

i.e. the field's own mean-plane downwash response divided by the analytic
horseshoe kernel's for the SAME sheet. That is the normalisation A13
reported as ~1.36x. No mu coupling, no Kutta row, no Newton drift: one
deterministic phi solve per grid.

  python tsd_a14_norm.py            default grid list
  python tsd_a14_norm.py 80 32 32   one grid
"""
import os
import sys

import numpy as np
from scipy.optimize import newton_krylov

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from wing3d.tsd import TSDSolver, flat_plate_jump     # noqa: E402

ALPHA = 1.25
MACH = float(os.environ.get('A14_NORM_M', '0.3'))
SPAN = 1.6


def phi_solve(s, epochs=2, tol=1e-9):
    nx, ny, nz = s.nx, s.ny, s.nz

    def unpack(y):
        phi = np.zeros((nx, ny, nz))
        phi[1:, :, :] = y.reshape(nx - 1, ny, nz)
        return phi

    def func(y):
        s.phi[:] = unpack(y)
        R = s.residual()
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        R[0, :, :] = 0
        return R[1:].ravel()

    y0 = s.phi[1:].ravel().copy()
    resn = float('nan')
    for _ in range(max(1, epochs)):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        y0 = newton_krylov(func, y0, iter=80, verbose=False,
                           f_tol=1e-6, f_rtol=tol)
        resn = float(np.abs(func(y0)).max())
    s.phi[:] = unpack(y0)
    return resn


def one(nx=80, ny=32, nz=32, epochs=2, wide=False):
    s = TSDSolver(MACH, alpha_deg=ALPHA, span=SPAN, nx=nx, ny=ny, nz=nz,
                  linear=False)
    s.farfield = None
    s.phi[:] = 0.0
    s.kutta_on = True
    s.kutta_free_mu = True
    s.kutta_bound = True
    s.kutta_thick_only = True
    s.kutta_kernel = True
    s.kutta_hybrid = True
    s.init_free_mu()
    nc, ng = s.n_ctrl, len(s.wing_jy)
    # exact flat-plate sheet, constant along the span (2D-like, the case
    # the kernel itself was validated on), pi*alpha*F(x)
    mu = (np.pi * np.radians(ALPHA) * flat_plate_jump(s.x_ctrl)[:, None]
          * np.ones((1, ng)))
    s.mu_b = mu.copy()
    s.wake_scale = 1.0
    resn = phi_solve(s, epochs=epochs)
    iu = np.asarray(s.wing_ix)
    jy = np.asarray(s.wing_jy)
    wu, wl, wav = s.w_plane(wide=wide)
    wf = np.asarray(wav[np.ix_(iu, jy)], float)
    L = s.dlm_obs_matrix()
    E = s.dlm_field_equivalent()
    we = (E @ mu.ravel()).reshape(len(iu), ng)
    wfc = (L @ wf.ravel()).reshape(nc, ng)
    j = ng // 2
    al = np.radians(ALPHA)
    wea, wfa = we[:, j] / al, wfc[:, j] / al
    m = np.abs(wea) > 0.25
    ratio = float(np.median(wfa[m] / wea[m])) if m.sum() else float('nan')
    zu, zl = s._stencil_heights()
    out = dict(nx=nx, ny=ny, nz=nz, nc=nc, ng=ng, x=s.x_ctrl.tolist(),
               wf=wfa.tolist(), we=wea.tolist(), ratio=ratio, res=resn,
               z_up=float(zu), z_lo=float(zl),
               dz_below=float(s.dz[s.k0 - 1]),
               dz_above=float(s.dz[min(s.k0, len(s.dz) - 1)]),
               z_up_over_dz=float(abs(zu) / max(s.dz[s.k0], 1e-30)),
               z_lo_over_dz=float(abs(zl) / max(s.dz[s.k0 - 1], 1e-30)))
    print('grid %3dx%3dx%3d  nc=%2d  z_up=%+.5f z_lo=%+.5f  dz(below)=%.4f '
          ' ratio=%.3f  |F|=%.1e' % (nx, ny, nz, nc, zu, zl, out['dz_below'],
                                     ratio, resn), flush=True)
    return out


def main(args):
    grids = [(40, 16, 16), (56, 32, 32), (80, 32, 24), (80, 32, 32),
             (80, 32, 48), (120, 32, 32), (160, 32, 32)] \
        if not args else [(int(args[0]), int(args[1]), int(args[2]))]
    res = [one(*g, epochs=int(os.environ.get('A14_NORM_EP', '2')))
           for g in grids]
    import json
    p = os.path.join(REPO, 'tests3d', 'a14_norm_M%.3f.json' % MACH)
    with open(p, 'w') as f:
        json.dump(res, f, indent=1)
    print('\nwrote', p)
    print('\n%-13s %4s %9s %9s %9s %7s' % ('grid', 'nc', 'z_up/c',
                                           'z_lo/c', 'dz_below', 'ratio'))
    for r in res:
        print('%-13s %4d %9.5f %9.5f %9.4f %7.3f' % (
            '%dx%dx%d' % (r['nx'], r['ny'], r['nz']), r['nc'], r['z_up'],
            r['z_lo'], r['dz_below'], r['ratio']))


if __name__ == '__main__':
    main(sys.argv[1:])
