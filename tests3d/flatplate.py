import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import solve, assemble, build_wake

# flat plate: zero-thickness wing, chord 1, span 6, modeled as TWO coincident
# surfaces (upper/lower quads) + edge closure, so thickness sources ~cancel
# and the test isolates doublet/wake/Kutta behavior.
nc, ns = 20, 12
xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, nc + 1)))  # cosine, LE->TE
ys = np.linspace(-3, 3, ns + 1)
verts = []
up_id, lo_id = {}, {}
for j, y in enumerate(ys):
    for i, x in enumerate(xs):
        up_id[(i, j)] = len(verts)
        verts.append([x, y, 0.0])
for j, y in enumerate(ys):
    for i, x in enumerate(xs):
        lo_id[(i, j)] = len(verts)
        verts.append([x, y, -1e-6])
faces = []
for j in range(ns):
    for i in range(nc):
        faces.append([up_id[(i, j)], up_id[(i + 1, j)],
                      up_id[(i + 1, j + 1)], up_id[(i, j + 1)]])
for j in range(ns):
    for i in range(nc):
        faces.append([lo_id[(i, j)], lo_id[(i, j + 1)],
                      lo_id[(i + 1, j + 1)], lo_id[(i + 1, j)]])
m = Mesh(np.array(verts), faces)
print('panels:', m.npanels)
# wake strips: TE pairs (upper TE panel i=nc-1, lower TE panel)
nup = (ns) * (nc)
info = []
for j in range(ns):
    up = j * nc + (nc - 1)
    lo = nup + j * nc + (nc - 1)
    pa = np.array([1.0, ys[j], 0.0])
    pb = np.array([1.0, ys[j + 1], 0.0])
    info.append((up, lo, pa, pb))
wakes = build_wake(m, info)
a = np.radians(4)
for km in ['doublet', 'potential']:
    A, Brow = assemble(m, wakes, kutta_mode=km)
    res = solve(m, wakes, [np.cos(a), 0, np.sin(a)], A=A, Brow=Brow,
                kutta_mode=km)
    F = (-res['cp'][:, None] * m.normal * m.area[:, None]).sum(axis=0)
    print(km, 'CL=%+.4f' % (F[2] / 6), '(thin-airfoil AR6 ~ +0.33)')
