"""A16 (non-planar geometry) shared harness.

One coupled Kutta solve on a low-RAM diagnostic grid, with hooks for
(a) a custom z-node distribution (plane-thickness / snapping study) and
(b) a future non-planar surface.  Uses tests3d/tsd_m_sup.run_sup so the
Kutta closures are exactly the ones the other agents measure.

Nothing here writes to cases/openfoam or the container.
"""
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from wing3d.tsd import TSDSolver, make_grid        # noqa: E402
from tsd_m_sup import run_sup                      # noqa: E402


def zgrid(nz, zlim=(-7.0, 7.0), core=(-0.5, 0.5), cluster=0.75):
    """The solver's own z distribution (make_grid's core(), excerpted).

    Same construction as wing3d.tsd.make_grid's `core(..., center_cluster)`
    so a default call reproduces the stock grid bit-for-bit; `cluster`
    and `core` are exposed so dz near the wing can be varied at fixed nz.
    """
    n = nz
    def core_fn(n, lo, hi, f0, f1, frac=0.6):
        nc = max(8, int(n * frac))
        u = np.linspace(0, 1, nc)
        A = cluster
        u = u + A * np.sin(2 * np.pi * u) / (2 * np.pi)
        c = f0 + (f1 - f0) * u
        nl = max(2, (n - nc) // 2)
        nr = max(2, n - nc - nl)
        left = lo + (f0 - lo) * np.geomspace(0.05, 1, nl) if nl > 1 else []
        right = f1 + (hi - f1) * np.geomspace(0.05, 1, nr) if nr > 1 else []
        g = np.unique(np.concatenate(
            [[lo], np.atleast_1d(left), c, np.atleast_1d(right), [hi]]))
        return np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(g)), g)
    return core_fn(nz, zlim[0], zlim[1], core[0], core[1])


def make_solver(nx=40, ny=16, nz=16, span=1.6, alpha=1.25, M=0.3,
                zdist=None, profile=None):
    """A bare TSDSolver with the stock grid, or `zdist` substituted."""
    s = TSDSolver(M, alpha_deg=alpha, span=span, nx=nx, ny=ny, nz=nz,
                  omega=0.8, profile=profile)
    if zdist is not None:
        z = np.asarray(zdist, float)
        assert z.shape == (nz,)
        s.z = z
        s.dz = np.diff(z)
        s.zc = 0.5 * (z[:-1] + z[1:])
        s.k0 = int(np.argmin(np.abs(z)))
    return s


def dzcore(s):
    """dz on the two faces bracketing the wing plane (fractions of chord)."""
    k = s.k0
    return (s.dz[k - 1], s.dz[min(k, len(s.dz) - 1)])


def run(s, M=0.3, epochs=1, alpha=1.25, span=None, profile=None,
        section=None, verbose=False):
    span = s.span if span is None else span
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=False,
                            profile=profile, section=section, base=s,
                            alpha_deg=alpha)
    return L


def span_loading(s):
    """gamma(y) per unit span from the mu sheet (mode 5/6/7) or gamma_b."""
    if getattr(s, 'kutta_free_mu', False) and s.mu_b.size:
        G = np.asarray(s.mu_circulation(), float)
    else:
        G = np.asarray(s.gamma_b, float)
    return np.asarray(s.yc[np.asarray(s.wing_jy)], float), G


def main(argv):
    nx, ny, nz = (int(v) for v in (argv[:3] if len(argv) >= 3
                                   else (40, 16, 16)))
    ep = int(os.environ.get('A16_EP', '1'))
    os.environ.setdefault('WING3D_KUTTA', '6')
    s = make_solver(nx, ny, nz)
    zu, zl = dzcore(s)
    L = run(s, epochs=ep)
    y, G = span_loading(s)
    print('g%dx%dx%d ep%d mode%s CL=%.5f CD=%.5f dz_below=%.4f '
          'dz_above=%.4f Gmid=%.5f' % (
              nx, ny, nz, ep, os.environ['WING3D_KUTTA'], L['CL'], L['CD'],
              zu, zl, G[len(G) // 2]))


if __name__ == '__main__':
    main(sys.argv[1:])
