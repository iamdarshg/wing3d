"""Inspect the actual LE point cloud for rae vs naca STL slices."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from wing3d.geometry import read_stl

for name in ('rae2822', 'naca0012'):
    mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\%s.stl'
                    % name)
    V = np.asarray(mesh.verts, float)
    m = np.abs(V[:, 1]) < 0.02
    P = V[m]
    near = P[P[:, 0] < 2e-3]
    print('== %s: slice pts=%d, pts with x<0.002: %d' % (name, len(P), len(near)))
    if len(near):
        order = np.argsort(near[:, 0])
        q = near[order]
        print('  x range of nose pts: %.6f .. %.6f' % (q[:, 0].min(), q[:, 0].max()))
        for xt in (0.0, 1e-4, 5e-4, 1e-3, 2e-3):
            sel = np.abs(q[:, 0] - xt) < 1e-4
            if sel.any():
                zz = np.sort(q[sel, 2])
                print('  x~%.4f: %d pts z=[%+.5f .. %+.5f]'
                      % (xt, sel.sum(), zz[0], zz[-1]))
    # absolute nose-tip: min-x vertex and its z
    imin = np.argmin(V[:, 0])
    print('  min-x vertex: (%.6f, y=%.4f, z=%+.6f)' %
          (V[imin, 0], V[imin, 1], V[imin, 2]))
