import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing, Mesh
from wing3d.primitives import ellipsoid, transform, merge
from wing3d.shapes import loft, rounded_rect
from wing3d.solver import assemble

scale = 0.2
L = 15.0 * scale
prof = [(0.00, 0.02, 0.0, 0.02), (0.06, 0.25, 0.0, 0.25),
        (0.15, 0.55, 0.05, 0.55), (0.30, 0.75, 0.1, 0.75),
        (0.45, 0.80, 0.1, 0.80), (0.60, 0.78, 0.1, 0.85),
        (0.75, 0.70, 0.1, 0.80), (0.90, 0.55, 0.1, 0.60),
        (1.00, 0.30, 0.1, 0.35)]
stations = []
for fx, hw, zc, hh in prof:
    loop = rounded_rect(hw * scale, hh * scale, 0.0, zc * scale, 16)
    stations.append((fx * L, loop))
fuse = loft(stations)
canopy = ellipsoid(a=2.2 * scale, b=0.45 * scale, c=0.45 * scale,
                   nlat=10, nlon=16, center=(0.32 * L, 0, 0.75 * scale))
wing = build_wing('0012', span=9.0 * scale, chord=3.2 * scale,
                  taper=0.25, sweep_deg=32.0, n_chord=20, n_span=10)
wing = transform(wing, translate=(0.45 * L, 0, -0.15 * scale))

for name, m in [('fuse', fuse), ('canopy', canopy), ('wing', wing)]:
    A, _ = assemble(m, [])
    print(name, 'npan=', m.npanels, 'cond=', round(float(np.linalg.cond(A)), 1))
m2 = merge([fuse, canopy])
A, _ = assemble(m2, [])
print('fuse+canopy npan=', m2.npanels, 'cond=', round(float(np.linalg.cond(A)), 1))
m3 = merge([fuse, wing])
A, _ = assemble(m3, [])
print('fuse+wing npan=', m3.npanels, 'cond=', round(float(np.linalg.cond(A)), 1))
