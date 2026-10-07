"""C-grid blockMeshDict generator for NACA0012 (2D transonic anchor)."""
import numpy as np

NACA = '0012'


def naca4(chord=1.0, n=121):
    m = p = 0.0
    t = 0.12
    beta = np.linspace(0, np.pi, n)
    x = (1 - np.cos(beta)) / 2
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                  + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    return x, yt


def emit(path, R=20.0, nx_surf=110, ny=36, nx_wake=50):
    x, yt = naca4(n=nx_surf + 1)
    te = 0.00126  # finite TE half-thickness (x=1)
    xu = x[::-1]
    yu = yt[::-1]
    xl = x
    yl = -yt
    with open(path, 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; '
                'object blockMeshDict; }\nconvertToMeters 1.0;\n')
        V = []
        # points: LE(0) upper TE(1) | farfield inlet-top(2) inlet-bot(3)
        # outlet-top(4) outlet-bot(5) | TE upper(6) TE lower(7)
        V.append((0, 0, 0))            # 0 LE
        V.append((1, te, 0))           # 1 TE upper
        V.append((1, -te, 0))          # 2 TE lower
        V.append((-R, R, 0))           # 3 far TL
        V.append((-R, -R, 0))          # 4 far BL
        V.append((1 + R, R, 0))        # 5 far TR
        V.append((1 + R, -R, 0))       # 6 far BR
        V.append((1 + R, te, 0))       # 7 wake top
        V.append((1 + R, -te, 0))      # 8 wake bot
        f.write('vertices\n(\n')
        for i, (px, py, pz) in enumerate(V):
            f.write(f'  ({px} {py} {pz}) // {i}\n')
        # z offset copies
        for i, (px, py, pz) in enumerate(V):
            f.write(f'  ({px} {py} 0.1) // {i + 9}\n')
        f.write(');\nblocks\n(\n')
        # upper passage: LE(0) TEu(1) farTR(5) farTL(3)
        f.write(f'  hex (0 1 5 3 9 10 14 12) ({nx_surf} {ny} 1) '
                f'simpleGrading (1 10 1)\n')
        # lower passage: LE(0) TEl(2) farBR(6) farBL(4)
        f.write(f'  hex (0 3 4 6 2 9 12 13 15 11) ({nx_surf} {ny} 1) '
                f'simpleGrading (1 0.1 1)\n')
        f.write(');\n')
        f.write('edges\n(\n')
        su = ' '.join(f'({a} {b} 0)' for a, b in zip(xu[::4], yu[::4]))
        sl = ' '.join(f'({a} {b} 0)' for a, b in zip(xl[::4], yl[::4]))
        f.write(f'  spline 0 1 ( {su} )\n')
        f.write(f'  spline 0 2 ( {sl} )\n')
        f.write(');\nboundary\n(\n')
        f.write('  airfoil { type wall; faces ( (0 1 10 9) (0 9 11 2) ); }\n')
        f.write('  farfield { type patch; faces ( (3 5 14 12) (4 6 15 13) '
                '(3 12 13 4) ); }\n')
        f.write('  outlet { type patch; faces ( (5 6 15 14) ); }\n')
        f.write('  frontback { type empty; faces ( (0 3 12 9) (1 5 14 10) '
                '(0 2 11 9) (2 6 15 11) ); }\n')
        f.write(');\nmergePatchPairs ( );\n')


if __name__ == '__main__':
    import sys
    emit(sys.argv[1] if len(sys.argv) > 1 else 'blockMeshDict')
    print('wrote C-grid dict')
