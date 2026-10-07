import numpy as np
import time
import sys
from wing3d.geometry import Mesh
from wing3d.solver import solve

nlat, nlon = int(sys.argv[1]), int(sys.argv[2])
verts = []
for i in range(nlat + 1):
    th = np.pi * i / nlat
    for j in range(nlon):
        ph = 2 * np.pi * j / nlon
        verts.append([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                      np.cos(th)])
verts = np.array(verts)
faces = []
for i in range(nlat):
    for j in range(nlon):
        a = i * nlon + j
        b = i * nlon + (j + 1) % nlon
        c = (i + 1) * nlon + j
        d = (i + 1) * nlon + (j + 1) % nlon
        if i > 0:
            faces.append([a, c, b])
        if i < nlat - 1:
            faces.append([b, c, d])
m = Mesh(verts, faces)
t0 = time.time()
res = solve(m, [], [1, 0, 0])
t1 = time.time()
cp = res['cp']
C = m.centroid
th = np.arccos(np.clip(C[:, 0], -1, 1))
exact = 1 - 2.25 * np.sin(th) ** 2
print(f'{nlat}x{nlon} npan={m.npanels} '
      f'rms={np.sqrt(((cp - exact) ** 2).mean()):.4f} '
      f'max={np.abs(cp - exact).max():.4f} t={t1 - t0:.1f}s')
