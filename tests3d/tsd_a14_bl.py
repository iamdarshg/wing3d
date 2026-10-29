"""A14 probe: the BL transpiration decomposition, no long solve.

Answers the two questions A13 left open on the viscous side without a full
ladder run:

  1. Does d(delta*)/dx reach `wav = (wu+wl)/2`, the only thing the mode-7
     row reads?  It does through the ASYMMETRIC part only
     (viscous camber (ddu-ddl)/2); the common source part (ddu+ddl)/2
     cancels, which is the CORRECT sign for a displacement thickness that
     thickens both surfaces.
  2. How big is each part, in alpha units, at mid span?

Usage:  python tsd_a14_bl.py [M] [epochs]
Runs one short coupled solve with WING3D_BL=1 on (needs a converged
state to march the BL on), then prints the decomposition.
"""
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_coupled import build_case                       # noqa: E402
from tsd_m_sup import run_sup                            # noqa: E402
from tsd_campaign import build_profile, build_section    # noqa: E402


def main(M=0.85, epochs=2, span=1.6, geo='dfvlr_r4'):
    os.environ.setdefault('WING3D_KUTTA', '7')
    os.environ.setdefault('WING3D_BL', '1')
    sec = build_section(geo)
    prof = build_profile(geo)
    base, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              span, section=sec, profile=prof)
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=False,
                            profile=prof, section=sec, base=base)
    bt = L['bl'].get('wav_terms')
    if not bt:
        print('no BL terms recorded')
        return
    al = np.radians(1.25)
    print('\nCL=%.5f   delta*/c max=%.3e' % (L['CL'], max(
        np.max(L['bl']['dstar_u']), np.max(L['bl']['dstar_l']))))
    print('x/c        :', ' '.join('%8.3f' % v for v in xc[::2]))
    print('ddu  (upper):', ' '.join('%8.2f' % v for v in
                                    np.asarray(bt['ddu_mid'])[::2]))
    print('ddl  (lower):', ' '.join('%8.2f' % v for v in
                                    np.asarray(bt['ddl_mid'])[::2]))
    print('source/2    :', ' '.join('%8.2f' % v for v in
                                    np.asarray(bt['half_mid'])[::2]))
    print('camber mean :', ' '.join('%8.2f' % v for v in
                                    np.asarray(bt['mean_mid'])[::2]))
    print('mean/alpha  :', ' '.join('%8.4f' % v for v in
                                    np.asarray(bt['wav_alpha_mid'])[::2]))


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.85,
         int(sys.argv[2]) if len(sys.argv) > 2 else 2)
