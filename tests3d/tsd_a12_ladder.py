"""A12 ladder: mode 6 = analytic horseshoe doublet-lattice tangency rows.

  python tsd_a12_ladder.py A     M0.3 symmetric   (test A)
  python tsd_a12_ladder.py Ac    M0.3 cambered dfvlr (test B, test C)
  python tsd_a12_ladder.py D     transonic dfvlr 0.85/0.925/1.00 (ENTROPY=1)
  python tsd_a12_ladder.py B     dfvlr M0.85 + WING3D_BL=1
  python tsd_a12_ladder.py F     free_k (mode 0) at M0.85

Same harness as tsd_a11_ladder.py (`one` imported from it) so the grids,
alpha, Re and the anchors are identical; only the closure flag moves.
"""
import os
import sys

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_a11_ladder import one    # noqa: E402  (shared harness)


def main(ladder):
    ep = int(os.environ.get('A12_EP', os.environ.get('A11_EP', '2')))
    if ladder == 'A':
        one(0.3, epochs=ep, use_profile=False, prefix='a12',
            note='mode6 M0.3 symmetric')
    elif ladder == 'Ac':
        one(0.3, epochs=ep, prefix='a12', note='mode6 M0.3 cambered')
    elif ladder == 'D':
        for M in (0.85, 0.925, 1.0):
            one(M, epochs=ep, prefix='a12', note='mode6 transonic')
    elif ladder == 'B':
        one(0.85, epochs=ep, prefix='a12', note='mode6 + BL')
    elif ladder == 'F':
        one(0.85, epochs=ep, free_k=True, prefix='a12', note='free_k legacy')
    else:
        raise SystemExit('unknown ladder %s' % ladder)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'A')
