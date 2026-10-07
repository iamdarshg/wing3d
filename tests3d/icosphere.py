import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import solve


def icosphere(ndiv=3):
    t = (1 + np.sqrt(5)) / 2
    v = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
                  [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
         (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
         (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
         (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    verts = [p for p in v]
    faces = [list(t) for t in f]
    for _ in range(ndiv):
        mid = {}

        def midpt(a, b):
            key = (min(a, b), max(a, b))
            if key not in mid:
                mid[key] = len(verts)
                verts.append((verts[a] + verts[b]) / 2)
            return mid[key]

        newf = []
        for (a, b, c) in faces:
            ab, bc, ca = midpt(a, b), midpt(b, c), midpt(c, a)
            newf += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = newf
    verts = np.array(verts)
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)
    return Mesh(verts, faces)


for ndiv in [2, 3]:
    m = icosphere(ndiv)
    res = solve(m, [], [1, 0, 0])
    cp = res['cp']
    C = m.centroid
    th = np.arccos(np.clip(C[:, 0], -1, 1))
    exact = 1 - 2.25 * np.sin(th) ** 2
    print(f'icosphere div={ndiv} npan={m.npanels} '
          f'rms={np.sqrt(((cp - exact) ** 2).mean()):.4f} '
          f'max={np.abs(cp - exact).max():.4f}')
