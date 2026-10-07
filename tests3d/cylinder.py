import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import solve

# circular cylinder extruded in y, span 8, radius 1, 64 around, 16 along
nr, ns = 64, 16
verts = []
for j in range(ns + 1):
    y = -4 + 8 * j / ns
    for i in range(nr):
        ph = 2 * np.pi * i / nr
        verts.append([np.cos(ph), y, np.sin(ph)])
verts = np.array(verts)
faces = []
for j in range(ns):
    for i in range(nr):
        i2 = (i + 1) % nr
        faces.append([j * nr + i, j * nr + i2,
                      (j + 1) * nr + i2, (j + 1) * nr + i])
# end caps (fans)
for j in (0, ns):
    c = np.array([0.0, -4.0 if j == 0 else 4.0, 0.0])
    ci = len(verts)
    verts = np.vstack([verts, c])
    ring = [j * nr + i for i in range(nr)]
    for i in range(nr):
        i2 = (i + 1) % nr
        faces.append([ring[i], ring[i2], ci])
m = Mesh(verts, faces)
print('panels:', m.npanels)
res = solve(m, [], [1, 0, 0])
cp = res['cp']
C = m.centroid
mid = np.abs(C[:, 1]) < 0.3
th = np.arctan2(C[mid][:, 2], C[mid][:, 0])  # angle from +x in x-z plane
exact = 1 - 4 * np.sin(th) ** 2
got = cp[mid]
print('midspan barrel: rms=%.4f max=%.4f' %
      (np.sqrt(((got - exact) ** 2).mean()), np.abs(got - exact).max()))
for deg in [0, 30, 60, 90]:
    k = int(np.argmin(abs(th * 180 / np.pi - deg)))
    print(f'  phi={deg}: num={got[k]:+.3f} exact={exact[k]:+.3f}')
