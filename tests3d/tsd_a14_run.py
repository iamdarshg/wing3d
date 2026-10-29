"""A14 run sets (low-RAM diagnostic grids).

The box is saturated, so every comparison here is RELATIVE on a reduced
grid; absolute CL is NOT comparable to the g80x32x32 anchors. Sets:

  python tsd_a14_run.py ref      M0.3 symmetric, ep4, mode 6 / mode 7 /
                                 mode 7 + B'                    (80x24x24)
  python tsd_a14_run.py ladder   dfvlr cambered, entropy=1, ep2,
                                 M0.85/0.925/1.00, B' off and on
  python tsd_a14_run.py grid     M0.3 symmetric, ep2, B' OFF -- the
                                 dz_core/nx scaling of the residual
  python tsd_a14_run.py bl       M0.3 symmetric, ep2, WING3D_BL=1
                                 with BL_SCALE 1 (physical) and 50
                                 (liveness/sign of dCL/d(delta*))
"""
import os
import sys

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_a14_ladder import one14          # noqa: E402

G = (80, 24, 24)
EP = int(os.environ.get('A14_EP', '4'))
EPL = int(os.environ.get('A14_EPL', '2'))


def ref():
    for km, tk in (('6', '0'), ('7', '0'), ('7', '1')):
        os.environ['WING3D_KUTTA'] = km
        os.environ['WING3D_KUTTA_THICK'] = tk
        os.environ['WING3D_THICK_EP'] = '1'
        one14(0.3, epochs=EP, use_profile=False, grid=G,
              note='ref set: mode %s thick-ref %s' % (km, tk))


def ladder():
    os.environ['WING3D_KUTTA'] = '7'
    os.environ['WING3D_ENTROPY'] = '1'
    os.environ['WING3D_THICK_EP'] = '2'
    for M in (0.85, 0.925, 1.0):
        for tk in ('0', '1'):
            os.environ['WING3D_KUTTA_THICK'] = tk
            one14(M, epochs=EPL, grid=(80, 32, 32),
                  note="ladder M%.3f thick-ref %s" % (M, tk))


def grid():
    os.environ['WING3D_KUTTA'] = '7'
    os.environ['WING3D_ENTROPY'] = '0'
    os.environ['WING3D_KUTTA_THICK'] = '0'
    for g in ((56, 32, 32), (80, 32, 24), (80, 32, 32), (80, 32, 48),
              (120, 32, 32)):
        one14(0.3, epochs=EPL, use_profile=False, grid=g,
              note='grid study %s' % (g,))


def bl():
    os.environ['WING3D_KUTTA'] = '7'
    os.environ['WING3D_ENTROPY'] = '0'
    os.environ['WING3D_KUTTA_THICK'] = '0'
    os.environ['WING3D_BL'] = '1'
    for sc in ('1', '50'):
        os.environ['WING3D_BL_SCALE'] = sc
        one14(0.3, epochs=EPL, use_profile=False, grid=G,
              note='BL scale %s' % sc)


if __name__ == '__main__':
    {'ref': ref, 'ladder': ladder, 'grid': grid, 'bl': bl}[
        sys.argv[1] if len(sys.argv) > 1 else 'ref']()
