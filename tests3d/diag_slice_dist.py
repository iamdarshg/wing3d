"""Inspect the STL slice point distribution (chordwise spacing)."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from wing3d.geometry import read_stl

mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
V = np.asarray(mesh.verts, float)
print('unique verts total:', len(V))
print('x extent %.4f .. %.4f  y %.4f .. %.4f  z %.4f .. %.4f'
      % (V[:, 0].min(), V[:, 0].max(), V[:, 1].min(), V[:, 1].max(),
         V[:, 2].min(), V[:, 2].max()))
# slice near mid-span
m = np.abs(V[:, 1]) < 0.02
P = V[m]
print('slice pts (|y|<0.02):', len(P))
ys = np.unique(np.round(P[:, 1], 4))
print('distinct y stations in slice:', len(ys), ys[:8], '...')
# chordwise distribution at one station
if len(ys):
    station = ys[len(ys) // 2]
    Q = P[np.abs(P[:, 1] - station) < 1e-4]
    xs = np.sort(Q[:, 0])
    print('station y=%.4f: %d pts, x %.4f..%.4f' % (station, len(Q), xs[0], xs[-1]))
    dx = np.diff(xs)
    print('dx min %.6f  max %.6f  median %.6f' % (dx.min(), dx.max(), np.median(dx)))
    for xq in [0.0005, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]:
        i = np.searchsorted(xs, xq)
        lo = max(0, i - 4)
        print('  x~%.4f  local dx = %s' % (xq, np.array2string(dx[lo:lo + 4], precision=5)))
