"""A17 diagnostic: does the stored phi carry a +-J/2 offset at the cut?

Runs ONE coupled Kutta-6 solve on the A16 diagnostic grid, then dumps the
z-structure of the field across the branch cut at mid chord / mid span so
the offset can be measured instead of assumed.

    python tmp_a17_diag.py [nx ny nz]
"""
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))
os.environ.setdefault('WING3D_KUTTA', '6')

from a16_common import make_solver, dzcore, run            # noqa: E402


def jump_at(s):
    """The bound-sheet jump the residual actually applies, (nix, njy)."""
    xc = np.clip(s.xc[s.wing_ix] / s.chord, 0.0, 1.0)
    if getattr(s, 'kutta_free_mu', False):
        jump = np.empty((len(s.wing_ix), len(s.wing_jy)))
        for p in range(len(s.wing_jy)):
            jump[:, p] = np.interp(xc, s.x_ctrl, s.mu_b[:, p])
        return jump
    if getattr(s, 'kutta_bound', False):
        from wing3d.tsd import flat_plate_jump
        return flat_plate_jump(xc)[:, None] * s.gamma_b[None, :]
    return np.asarray(s.wing_jump)


def main(argv):
    nx, ny, nz = (int(v) for v in (argv if len(argv) >= 3 else (40, 24, 16)))
    s = make_solver(nx, ny, nz)
    L = run(s, epochs=int(os.environ.get('A17_EP', '1')))
    k0 = s.k0
    print('\nCL=%.5f  nx,ny,nz=%d,%d,%d  k0=%d z[k0]=%.6f'
          % (L['CL'], s.nx, s.ny, s.nz, k0, s.z[k0]))
    print('z[k0-3..k0+3] =', np.round(s.z[k0 - 3:k0 + 4], 5))
    print('dz  around k0 =', np.round(s.dz[k0 - 3:k0 + 3], 5))

    J = jump_at(s)
    iu = np.asarray(s.wing_ix)
    jy = np.asarray(s.wing_jy)
    jm = len(jy) // 2
    j = jy[jm]
    print('\n--- column j=%d (mid span), i = wing chord nodes ---' % j)
    print('%6s %9s %9s %9s %9s %9s' % ('xc', 'J', 'phi[k0]-1',
                                       'dv/dz', 'dvz-J/dz', 'dphi/dz above'))
    for n, i in enumerate(iu):
        fc = s.phi[i, j, k0] - s.phi[i, j, k0 - 1]
        dzf = s.z[k0] - s.z[k0 - 1]
        up = (s.phi[i, j, k0 + 1] - s.phi[i, j, k0]) / (s.z[k0 + 1] - s.z[k0])
        print('%6.3f %9.5f %9.5f %9.5f %9.5f %9.5f'
              % (s.xc[i] / s.chord, J[n, jm], fc, fc / dzf,
                 fc / dzf - J[n, jm] / dzf, up))

    # lateral: does the offset die outside the wing span?
    print('\n--- row i=mid chord, across span ---')
    im = iu[len(iu) // 2]
    n0 = len(iu) // 2
    for m, jj in enumerate(jy):
        fc = s.phi[im, jj, k0] - s.phi[im, jj, k0 - 1]
        print('  y=%7.3f  J=%9.5f  phi[k0]-phi[k0-1]=%9.5f  fc/dz=%9.5f'
              % (s.yc[jj], J[n0, m], fc, fc / (s.z[k0] - s.z[k0 - 1])))

    # what the estimators actually read at mid chord/mid span
    kp, km = min(k0 + 1, s.nz - 1), max(k0 - 2, 0)
    print('\nwgx by plane at (xc=0.75, y=mid):')
    i4 = int(iu[int(0.75 * (len(iu) - 1))])
    for k in range(max(0, k0 - 3), min(s.nz, k0 + 4)):
        g = (s.phi[i4 + 1, j, k] - s.phi[i4 - 1, j, k]) / (
            s.x[i4 + 1] - s.x[i4 - 1])
        print('  k=%2d z=%+.4f  dphi/dx=%+.6f  (gradx+1=%+.6f)'
              % (k, s.z[k], g, 1 + g))
    wu, wl, wav = s.w_plane()
    print('w_plane at (xc=0.75,mid): wu=%+.5f wl=%+.5f wav=%+.5f (alpha=%.4f)'
          % (wu[i4, j], wl[i4, j], wav[i4, j], s.alpha))
    out = {'nz': nz, 'CL': L['CL'], 'k0': k0, 'z0': float(s.z[k0]),
           'dz_below': float(s.dz[k0 - 1]), 'dz_above': float(s.dz[k0])}
    import json
    with open(os.path.join(REPO, 'tmp_a17_diag_%d.json' % nz), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main(sys.argv[1:])
