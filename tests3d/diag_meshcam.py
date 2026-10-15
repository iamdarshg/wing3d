"""Extract the wing section mean line from an OpenFOAM polyMesh points file.

Filters the y=0 slab, bins x, takes (zmax+zmin)/2 = mean line; reports the
camber amplitude, so we can tell WHICH STL a mesh was built from.
"""
import sys
import numpy as np


def load_points(path):
    lines = open(path, 'r').read().splitlines()
    i = 0
    for i, l in enumerate(lines):
        if l.strip() == '(':
            break
    vals = []
    for l in lines[i + 1:]:
        t = l.strip()
        if not t or t == ')':
            break
        p = t.replace('(', ' ').replace(')', ' ').split()
        vals.append((float(p[0]), float(p[1]), float(p[2])))
    return np.array(vals)


def camber_of(path, tag):
    P = load_points(path)
    sel = (np.abs(P[:, 1]) < 0.05) & (np.abs(P[:, 2]) < 0.2) & \
          (P[:, 0] > -0.05) & (P[:, 0] < 1.05)
    S = P[sel]
    xs = np.linspace(0.0, 1.0, 41)
    mean = np.full(len(xs), np.nan)
    thick = np.full(len(xs), np.nan)
    for k in range(len(xs) - 1):
        m = (S[:, 0] >= xs[k]) & (S[:, 0] < xs[k + 1])
        if m.sum() < 4:
            continue
        z = S[m, 2]
        up = z[z > 0.0]
        lo = z[z < 0.0]
        if len(up) == 0 or len(lo) == 0:
            continue
        zu = up.min()   # first mesh point above the upper surface
        zl = lo.max()   # first mesh point below the lower surface
        mean[k] = 0.5 * (zu + zl)
        thick[k] = zu - zl
    good = np.isfinite(mean)
    if good.sum() < 5:
        print('%-22s NO SECTION FOUND (%d pts in slab)' % (tag, len(S)))
        return
    m0 = mean[good] - mean[good].mean()
    print('%-22s nPts=%7d  camber_amp=%.4f  thick_max=%.4f  '
          'mean-line: %s' % (tag, len(S), float(np.max(np.abs(m0))),
                             float(np.nanmax(thick)),
                             ' '.join('%.3f' % v for v in
                                      m0[::max(1, good.sum() // 8)])))


for p, tag in [
        (sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'mesh')]:
    camber_of(p, tag)
