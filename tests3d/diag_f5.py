import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d import shapes as S
from wing3d.primitives import ellipsoid, wedge
from collections import Counter

scale = 0.1
L = 14.45 * scale
prof = [(0.00, 0.03, 0.0, 0.03), (0.05, 0.28, 0.0, 0.30),
        (0.12, 0.50, 0.02, 0.55), (0.25, 0.72, 0.05, 0.80),
        (0.40, 0.85, 0.05, 0.95), (0.60, 0.90, 0.05, 1.00),
        (0.75, 0.88, 0.05, 0.95), (0.88, 0.75, 0.08, 0.75),
        (0.96, 0.55, 0.10, 0.50), (1.00, 0.35, 0.10, 0.32)]
R = 1.0 * scale
stations = []
for fx, hw, zc, hh in prof:
    stations.append((fx * L, S.rounded_rect(hw * R, hh * R, 0.0, zc * R,
                                            16)))
fuse = S.loft(stations)
wing = S._place_wing(dict(code='0012', span=8.1 * scale, chord=2.9 * scale,
                          taper=0.30, sweep_deg=24.0, n_chord=14,
                          n_span=10, cosine=False),
                     (0.42 * L, 0, -0.35 * scale))
lex_r = wedge(x0=0.30 * L, x1=0.52 * L, y0=0.10 * scale,
              y1=1.15 * scale, z0=-0.30 * scale, z1=0.05 * scale,
              nx=6, ny=4)


def open_edges(mesh):
    edges = Counter()
    for f in mesh.faces:
        m = len(f)
        V = [tuple(np.round(mesh.verts[k], 9)) for k in f]
        for k in range(m):
            edges[tuple(sorted([V[k], V[(k + 1) % m]]))] += 1
    return sum(1 for c in edges.values() if c == 1)


for name, m in [('fuse', fuse), ('wing', wing), ('lex_r', lex_r)]:
    print('%s: panels=%d open_edges=%d' % (name, m.npanels, open_edges(m)))
