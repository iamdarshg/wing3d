"""A16: which piece of the non-planar path is usable at 5 deg dihedral?

Splits the implementation into its three independent mechanisms and runs
each combination: the sheet/cut face (WING3D_NP_CUT), the wake face
(WING3D_NP_WAKE), both vs the planar reference. Reported at a dihedral
angle whose VLM effect is -0.31% (so a usable path must stay within a few
percent of the flat case, not collapse).
"""
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from a16_common import make_solver, run, span_loading      # noqa
import a16_nonplanar as NP                                 # noqa

GRID = (40, 24, 16)
DEG = float(os.environ.get('A16_DEG', '5'))


def one(tag, zsurf, cut='0', wake='0'):
    os.environ['WING3D_NP_CUT'] = cut
    os.environ['WING3D_NP_WAKE'] = wake
    s = make_solver(*GRID)
    if zsurf is not None:
        s.set_surface(zsurf)
    L = run(s, epochs=1)
    y, G = span_loading(s)
    clc = 2.0 * float(np.sum(np.asarray(G, float)
                             * s.dy[np.asarray(s.wing_jy)])) / s.S_ref()
    print('%-28s CL=%+.5f  CL_circ=%+.5f' % (tag, L['CL'], clc), flush=True)
    return L['CL'], clc


if __name__ == '__main__':
    os.environ['WING3D_KUTTA'] = os.environ.get('A16_MODE', '6')
    s = make_solver(*GRID)
    z = NP._wing_zsurf((s.x, s.y), 'dih', deg=DEG)
    print('--- flat reference ---')
    flat, flatc = one('flat', None)
    print('--- dihedral %.1f deg (cut x wake) ---' % DEG)
    for cut in ('0', '1'):
        for wake in ('0', '1'):
            one('dih cut=%s wake=%s' % (cut, wake), z, cut, wake)
