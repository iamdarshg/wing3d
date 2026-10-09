import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d import shapes as S
from wing3d.primitives import ellipsoid, wedge, merge
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces

scale = 0.1
L = 14.45 * scale
prof = [(0.00, 0.03, 0.0, 0.03), (0.05, 0.28, 0.0, 0.30),
        (0.12, 0.50, 0.02, 0.55), (0.25, 0.72, 0.05, 0.80),
        (0.40, 0.85, 0.05, 0.95), (0.60, 0.90, 0.05, 1.00),
        (0.75, 0.88, 0.05, 0.95), (0.88, 0.75, 0.08, 0.75),
        (0.96, 0.55, 0.10, 0.50), (1.00, 0.35, 0.10, 0.32)]
R = 1.0 * scale
stations = [(fx * L, S.rounded_rect(hw * R, hh * R, 0.0, zc * R, 16))
            for fx, hw, zc, hh in prof]
fuse = S.loft(stations)
canopy = ellipsoid(a=1.6 * scale, b=0.35 * scale, c=0.45 * scale,
                   nlat=10, nlon=16, center=(0.30 * L, 0, 0.55 * scale))
wing = S._place_wing(dict(code='0012', span=8.1 * scale, chord=2.9 * scale,
                          taper=0.30, sweep_deg=24.0, n_chord=14,
                          n_span=10, cosine=False),
                     (0.42 * L, 0, -0.35 * scale))
lex_r = wedge(x0=0.30 * L, x1=0.52 * L, y0=0.10 * scale,
              y1=1.15 * scale, z0=-0.30 * scale, z1=0.05 * scale,
              nx=6, ny=4)
V = np.array([1.0, 0, 0.0])
sref = 17.28 * scale * scale
base = [fuse, canopy, wing]
for name, extra in [('+lex_r', [lex_r])]:
    m = merge(base + extra)
    res = solve(m, [], V)
    f = pressure_forces(m, res['cp'], V, sref=sref)
    print('%s: npan=%d CL=%+.4f CDp=%.5f maxCp=%.1f' % (
        name, m.npanels, f['CL'], f['CDp'], np.abs(res['cp']).max()))
# wing + wakes (no lex/tails); wing offset = 176 + 288 = 464
ww, skipped = S._surface_wakes(wing, 464, cut=('y', 0.95 * scale))
m = merge(base)
res = solve(m, ww, V)
f = pressure_forces(m, res['cp'], V, sref=sref)
print('+wingwakes: CL=%+.4f CDp=%.5f maxCp=%.1f' % (
    f['CL'], f['CDp'], np.abs(res['cp']).max()))
