"""Launch checklist: verify case inflow matches intent (M, alpha).

Reads a LOCAL case dir (0/U + system/controlDict); asserts:
|U| == M*a_inf, angle(U) == alpha, magUInf == |U|.
Prevents mislabeled runs (2026-10-08 audit). Usage:
  python tests3d/of_verify.py cases/openfoam/transonic 0.8 1.25
"""
import re
import sys


def parse_field(path, key):
    txt = open(path).read()
    m = re.search(r'internalField\s+uniform\s+\(([^)]+)\)', txt)
    if not m:
        raise ValueError('no internalField in %s' % path)
    return [float(x) for x in m.group(1).split()]


def parse_mag(path):
    txt = open(path).read()
    m = re.search(r'magUInf\s+([\d.eE+-]+)', txt)
    return float(m.group(1)) if m else None


if __name__ == '__main__':
    root, M, alpha = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    import numpy as np
    U = parse_field(root + '/0/U', 'U')
    mag = parse_mag(root + '/system/controlDict')
    a_inf = np.sqrt(1.4 * 287.0 * 300.0)
    Uexp = M * a_inf
    ang = np.degrees(np.arctan2(U[2], U[0]))
    ok = True
    for name, got, exp, tol in [('|U|', np.linalg.norm(U), Uexp, 0.5),
                                ('alpha', ang, alpha, 0.05),
                                ('magUInf', mag, Uexp, 0.5)]:
        good = abs(got - exp) <= tol
        ok &= good
        print('%s: got %.3f exp %.3f %s' % (name, got, exp,
                                            'OK' if good else 'MISMATCH'))
    sys.exit(0 if ok else 1)
