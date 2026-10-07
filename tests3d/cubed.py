import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import solve


def cubed_sphere(n):
    verts = []
    faces = []
    # 6 faces of cube, gnomic projection to sphere
    for f in range(6):
        base = len(verts)
        us = np.linspace(-1, 1, n + 1)
        idx = {}
        for i, u in enumerate(us):
            for j, v in enumerate(us):
                if f == 0:
                    p = np.array([1., u, v])
                elif f == 1:
                    p = np.array([-1., u, v])
                elif f == 2:
                    p = np.array([u, 1., v])
                elif f == 3:
                    p = np.array([u, -1., v])
                elif f == 4:
                    p = np.array([u, v, 1.])
                else:
                    p = np.array([u, v, -1.])
                p = p / np.linalg.norm(p)
                idx[(i, j)] = len(verts)
                verts.append(p)
        for i in range(n):
            for j in range(n):
                faces.append([idx[(i, j)], idx[(i + 1, j)],
                              idx[(i + 1, j + 1)], idx[(i, j + 1)]])
    return Mesh(np.array(verts), faces)


for n in [6, 10]:
    m = cubed_sphere(n)
    res = solve(m, [], [1, 0, 0])
    cp = res['cp']
    C = m.centroid
    th = np.arccos(np.clip(C[:, 0] / np.linalg.norm(C, axis=1), -1, 1))
    exact = 1 - 2.25 * np.sin(th) ** 2
    print(f'cubed n={n} npan={m.npanels} '
          f'rms={np.sqrt(((cp - exact) ** 2).mean()):.4f} '
          f'max={np.abs(cp - exact).max():.4f}')
