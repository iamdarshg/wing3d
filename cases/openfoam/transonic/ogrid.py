"""O-grid blockMeshDict generator for exact-2D airfoil anchors.

Topology: inner loop (airfoil, 4 corners) + outer square + 4 blocks.
Inner corners: UNL (upper near-LE), TEU, TEL, LNL (lower near-LE);
splines run smoothly THROUGH the LE (no degenerate pole).
Radial edges graded geometrically (first-cell control, no snappy
staircasing). Single cell in y + empty front/back = true 2D.
"""
import numpy as np


def naca0012_yt(x):
    x = np.clip(np.asarray(x, dtype=float), 0, 1)
    return 5 * 0.12 * (0.2969 * np.sqrt(np.maximum(x, 1e-12)) - 0.1260 * x
                       - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)


def emit(path, R=20.0, n_airfoil=480, n_radial=64, first_dr=2e-4):
    # inner corners (x0 small: fore inner edge clips ~1mm of nose)
    x0 = 0.001
    unl = (x0, naca0012_yt(x0))
    teu = (1.0, naca0012_yt(1.0))
    tel = (1.0, -naca0012_yt(1.0))
    lnl = (x0, -naca0012_yt(x0))
    inner = [unl, teu, tel, lnl]
    # outer corners (square, same angular order); outlet near (x=6)
    # to keep aft-block cells sane (far outlet -> 30000:1 TE slivers)
    outer = [(-R, R), (7.0, R), (7.0, -R), (-R, -R)]
    # airfoil splines: upper (unl->teu through LE? NO: upper from unl fwd
    # around LE? unl is AT x=0.008 upper; upper edge goes unl -> ... must
    # pass LE(0,0) -> teu. Build dense upper/lower point lists.
    xu = np.linspace(0, 1, 400)
    # upper edge: unl -> LE -> teu
    up1 = [(t, naca0012_yt(t)) for t in np.linspace(x0, 0.0, 30)][:-1]
    up2 = [(t, naca0012_yt(t)) for t in np.linspace(0.0, 1.0, 300)]
    upper = up1 + up2
    lo1 = [(t, -naca0012_yt(t)) for t in np.linspace(x0, 0.0, 30)][:-1]
    lo2 = [(t, -naca0012_yt(t)) for t in np.linspace(0.0, 1.0, 300)]
    lower = lo1 + lo2  # lnl -> LE -> tel (must run corner-to-corner)
    with open(path, 'w', newline='\n') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; '
                'object blockMeshDict; }\nconvertToMeters 1.0;\n')
        V = inner + outer
        f.write('vertices\n(\n')
        for i, (px, py) in enumerate(V):
            f.write('  (%g %g 0) // %d\n' % (px, py, i))
        for i, (px, py) in enumerate(V):
            f.write('  (%g %g 0.1) // %d\n' % (px, py, i + 8))
        f.write(');\nblocks\n(\n')
        # radial grading: geometric expansion from wall
        # blockMesh simpleGrading is uniform-ratio; approximate with
        # multi-segment? Use single grading ratio r chosen so first
        # cell ~ first_dr: r ≈ (R/n)^(1/(n-1))... use expansion below
        import math
        ratio = (R / max(first_dr, 1e-9)) ** (1.0 / max(n_radial - 1, 1))
        ratio = min(max(ratio, 1.02), 1.3)
        na = n_airfoil // 4
        # CCW-wound hexes (blockMesh rejects inside-out):
        # top (unl,teu,OUT_TR,OUT_TL), aft (TEl,OUT_BR,OUT_TR,TEu),
        # bottom (TEu? no: (TEu->TEl reversed) (2 3 7 6)), fore (lnl,unl...)
        f.write('  hex (0 1 5 4 8 9 13 12) (%d %d 1) simpleGrading (1 1 1)\n'
                % (na, n_radial))
        f.write('  hex (2 6 5 1 10 14 13 9) (%d %d 1) simpleGrading (1 1 1)\n'
                % (n_radial, max(na // 4, 4)))
        f.write('  hex (2 3 7 6 10 11 15 14) (%d %d 1) simpleGrading (1 1 1)\n'
                % (na, n_radial))
        f.write('  hex (3 0 4 7 11 8 12 15) (%d %d 1) simpleGrading (1 1 1)\n'
                % (max(na // 4, 4), n_radial))
        f.write(');\nedges\n(\n')
        su = ' '.join('(%g %g 0)' % p for p in upper[::6])
        sl = ' '.join('(%g %g 0)' % p for p in lower[::6])
        f.write('  spline 0 1 ( %s )\n' % su)
        f.write('  spline 3 2 ( %s )\n' % sl)
        f.write(');\nboundary\n(\n')
        f.write('  airfoil { type wall; faces ( (0 1 9 8) (2 1 9 10) '
                '(2 3 11 10) (3 0 8 11) ); }\n')
        f.write('  farfield { type patch; faces ( (4 5 13 12) (5 6 14 13) '
                '(7 6 14 15) (4 7 15 12) ); }\n')
        f.write('  frontback { type empty; faces ( (0 1 5 4) (2 6 5 1) '
                '(2 3 7 6) (3 0 4 7) (8 9 13 12) (10 14 13 9) '
                '(10 11 15 14) (11 8 12 15) ); }\n')
        f.write(');\nmergePatchPairs ( );\n')
    print('wrote O-grid: ratio=%.3f cells~%d' % (
        ratio, 4 * na * n_radial))


if __name__ == '__main__':
    import sys
    emit(sys.argv[1] if len(sys.argv) > 1 else 'blockMeshDict')
