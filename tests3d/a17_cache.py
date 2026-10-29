"""A17: solve each case ONCE and cache the converged state (dict .npz).

Re-running solves per estimator costs 45 s each; the estimators are pure
post-processing of (phi, z, mu_b, gamma_b), so cache them once.

    python tmp_a17_solve.py          # writes tmp_a17_case_<tag>.npz
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))
os.environ.setdefault('WING3D_KUTTA', '6')

from a16_common import make_solver, run, zgrid      # noqa: E402

CASES = [('nz12', dict(nz=12)), ('nz16', dict(nz=16)),
         ('nz24', dict(nz=24)), ('nz32', dict(nz=32))]


def main():
    z0 = zgrid(16)
    for frac in (0.25, 0.5):
        CASES.append(('zoff%d' % int(frac * 100),
                      dict(nz=16, zdist=z0 + frac * 0.0406)))
    out = {}
    for tag, kw in CASES:
        f = os.path.join(REPO, 'tmp_a17_case_%s.npz' % tag)
        if os.path.exists(f) and os.environ.get('A17_FORCE') != '1':
            print('cached', tag, flush=True)
            continue
        s = make_solver(40, 24, kw['nz'], zdist=kw.get('zdist'))
        L = run(s, epochs=int(os.environ.get('A17_EP', '1')), verbose=False)
        np.savez_compressed(
            f, phi=s.phi, x=s.x, y=s.y, z=s.z, dx=s.dx, dy=s.dy, dz=s.dz,
            k0=s.k0, wing_ix=np.asarray(s.wing_ix),
            wing_jy=np.asarray(s.wing_jy), mu_b=s.mu_b,
            gamma_b=np.asarray(s.gamma_b), x_ctrl=s.x_ctrl,
            xc=s.xc, yc=s.yc, thick=s.thick, Minf=s.Minf, gamma=s.gamma,
            alpha=s.alpha, chord=s.chord, span=s.span,
            free_mu=np.array([bool(getattr(s, 'kutta_free_mu', False))]),
            CL_fp=np.array([L['CL']]), tag=np.array([tag]))
        nmid = len(s.wing_jy) // 2
        G = (s.mu_circulation() if getattr(s, 'kutta_free_mu', False)
             else np.asarray(s.gamma_b))
        out[tag] = {'CL_fp': float(L['CL']),
                    'CL_circ': 2.0 * float(np.sum(
                        G * s.dy[np.asarray(s.wing_jy)])) / (s.chord * s.span),
                    'z0': float(s.z[s.k0]),
                    'dz': [float(s.dz[s.k0 - 1]), float(s.dz[s.k0])],
                    'G_mid': float(G[nmid])}
        print('%-6s CL_fp=%.5f CL_circ=%.5f z0=%+.4f dz=%.4f/%.4f'
              % (tag, out[tag]['CL_fp'], out[tag]['CL_circ'],
                 out[tag]['z0'], *out[tag]['dz']), flush=True)
        del s
    with open(os.path.join(REPO, 'tmp_a17_cases.json'), 'w') as fh:
        json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
