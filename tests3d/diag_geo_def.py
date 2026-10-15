"""Definitive geometry reconciliation: campaign STLs vs masters vs tables vs meshes.

For each source prints the mean line (camber) at fixed stations, max
camber and its location, and max thickness. Surface-only sources (STL)
use zmax/zmin per bin; volume sources (mesh) use the first-crossing
method (nearest point above/below z=0).
"""
import hashlib
import struct
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

STL_DIR = r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl'
CAMP = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'


def read_stl(path):
    b = open(path, 'rb').read()
    n = struct.unpack('<I', b[80:84])[0]
    if 84 + 50 * n == len(b):
        f = np.frombuffer(b, dtype=np.uint8, count=50 * n, offset=84)
        f = f.reshape(n, 50)[:, :48].copy().view('<f4').reshape(n, 12)
        return f[:, 3:12].reshape(-1, 3)
    # ascii
    v = []
    for l in b.decode('utf-8', 'ignore').splitlines():
        t = l.strip().split()
        if t and t[0] == 'vertex':
            v.append([float(t[1]), float(t[2]), float(t[3])])
    return np.array(v)


def h5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()[:12]


def camber_surface(P, tag):
    """STL surface point cloud: zmax/zmin per x-bin = true surfaces."""
    sel = P[(np.abs(P[:, 1]) < 0.06) & (np.abs(P[:, 2]) < 0.4)]
    if len(sel) < 50:
        print('%-34s NO SECTION (%d pts)' % (tag, len(sel)))
        return
    return _report(sel[:, 0], sel[:, 2], tag)


def camber_mesh(P, tag):
    """Volume mesh points: first point above/below z=0 = surface."""
    sel = P[(np.abs(P[:, 1]) < 0.05) & (np.abs(P[:, 2]) < 0.4) &
            (P[:, 0] > -0.05) & (P[:, 0] < 1.05)]
    if len(sel) < 50:
        print('%-34s NO SECTION (%d pts)' % (tag, len(sel)))
        return
    xs = np.linspace(0, 1, 41)
    cm = np.full(40, np.nan)
    th = np.full(40, np.nan)
    for k in range(40):
        m = (sel[:, 0] >= xs[k]) & (sel[:, 0] < xs[k + 1])
        if m.sum() < 2:
            continue
        z = sel[m, 2]
        up = z[z > 1e-6]
        lo = z[z < -1e-6]
        if not len(up) or not len(lo):
            continue
        cm[k] = 0.5 * (up.min() + lo.max())
        th[k] = up.min() - lo.max()
    _print(cm, th, xs, tag)


def _report(x, z, tag):
    xs = np.linspace(0, 1, 41)
    cm = np.full(40, np.nan)
    th = np.full(40, np.nan)
    for k in range(40):
        m = (x >= xs[k]) & (x < xs[k + 1])
        if m.sum() < 2:
            continue
        zz = z[m]
        cm[k] = 0.5 * (zz.max() + zz.min())
        th[k] = zz.max() - zz.min()
    _print(cm, th, xs, tag)


def _print(cm, th, xs, tag):
    xs = np.linspace(0, 1, len(cm) + 1)
    xc = 0.5 * (xs[:-1] + xs[1:])
    good = np.isfinite(cm)
    if good.sum() < 5:
        print('%-34s unreadable' % tag)
        return
    st = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    vals = [np.interp(s, xc[good], cm[good]) for s in st]
    imax = np.nanargmax(np.abs(cm))
    print('%-34s cmax=%+.4f@x=%.2f | tmax=%.4f | camber@[.2 .3 .4 .5 .6 .7 .8 .9]=' %
          (tag, cm[imax], xc[imax], np.nanmax(th)) +
          ' '.join('%+.3f' % v for v in vals))


def load_of_points(path):
    lines = open(path).read().splitlines()
    i = 0
    for i, l in enumerate(lines):
        if l.strip() == '(':
            break
    v = []
    for l in lines[i + 1:]:
        t = l.strip()
        if not t or t == ')':
            break
        p = t.replace('(', ' ').replace(')', ' ').split()
        v.append((float(p[0]), float(p[1]), float(p[2])))
    return np.array(v)


def table_section(fn, args, tag):
    out = fn(*args)
    if len(out) == 2:
        xg, yg = out
        _report(np.hstack([xg, xg]), np.hstack([yg, -yg]), tag)
        return
    xu, yu, xl, yl = [np.asarray(a, float) for a in out]
    if xu[0] > xu[-1]:
        xu, yu = xu[::-1], yu[::-1]
    if xl[0] > xl[-1]:
        xl, yl = xl[::-1], yl[::-1]
    xg = np.linspace(0, 1, 400)
    _report(xg, np.interp(xg, xu, yu), tag + ' [upper]')
    _report(xg, np.interp(xg, xl, yl), tag + ' [lower]')
    cm = 0.5 * (np.interp(xg, xu, yu) + np.interp(xg, xl, yl))
    th = np.interp(xg, xu, yu) - np.interp(xg, xl, yl)
    _print(cm, th, xg[:400], tag + ' [meanline]')


def main():
    print('== hashes (md5[:12])')
    for p, tag in [
            (STL_DIR + r'\rae2822.stl', 'master rae'),
            (r'D:\CodeProjects\cfd\tmp_rae_case.stl', 'case rae (M0.850)'),
            (STL_DIR + r'\dfvlr_r4.stl', 'master dfvlr'),
            (r'D:\CodeProjects\cfd\tmp_dfvlr_case.stl', 'case dfvlr (M0.850)'),
            (STL_DIR + r'\naca0012.stl', 'master naca'),
            (CAMP + r'\naca0012_M0.850_a1.25\constant\triSurface\naca0012.stl',
             'case naca'),
            (STL_DIR + r'\biconvex.stl', 'master biconvex'),
            (CAMP + r'\biconvex_M0.850_a1.25\constant\triSurface\biconvex.stl',
             'case biconvex')]:
        try:
            print('  %-24s %s' % (tag, h5(p)))
        except OSError as e:
            print('  %-24s MISSING (%s)' % (tag, e))

    print('\n== STL sections')
    for p, tag in [
            (STL_DIR + r'\rae2822.stl', 'master rae stl'),
            (r'D:\CodeProjects\cfd\tmp_rae_case.stl', 'case rae stl'),
            (STL_DIR + r'\dfvlr_r4.stl', 'master dfvlr stl'),
            (r'D:\CodeProjects\cfd\tmp_dfvlr_case.stl', 'case dfvlr stl')]:
        try:
            P = read_stl(p)
            g = P[np.abs(P[:, 1]) < 0.06]
            print('  [%s] verts=%d yrange=(%.3f,%.3f) zrange=(%.3f,%.3f)'
                  % (tag, len(P), P[:, 1].min(), P[:, 1].max(),
                     P[:, 2].min(), P[:, 2].max()))
            camber_surface(P, '  ' + tag)
        except OSError:
            print('  [%s] MISSING' % tag)

    print('\n== coordinate tables')
    from gen_transonic_stl import (rae2822_xy, dfvlr_r4_xy, naca0012_xy,
                                   biconvex_xy)
    table_section(rae2822_xy, (), 'table rae')
    table_section(dfvlr_r4_xy, (), 'table dfvlr')
    xa, yt = naca0012_xy(n=400)
    _report(xa, yt, 'table naca [upper]')
    xb, yb = biconvex_xy(n=400)
    _report(xb, yb, 'table biconvex [upper]')

    print('\n== meshes')
    for p, tag in [(r'D:\CodeProjects\cfd\tmp_rae_new_pts', 'mesh rae'),
                   (r'D:\CodeProjects\cfd\tmp_dfvlr_pts', 'mesh dfvlr'),
                   (r'D:\CodeProjects\cfd\tmp_naca_pts', 'mesh naca'),
                   (r'D:\CodeProjects\cfd\tmp_bic_pts', 'mesh biconvex')]:
        try:
            camber_mesh(load_of_points(p), tag)
        except OSError:
            print('%-34s MISSING' % tag)


if __name__ == '__main__':
    main()
