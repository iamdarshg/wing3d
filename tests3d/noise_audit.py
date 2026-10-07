"""Noise floor + weak-spot audit: exact residual levels per geometry/condition."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from wing3d.geometry import Mesh, build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces
from wing3d.primitives import box, ellipsoid
from analytic_sphere import sphere_mesh

print('=== A. symmetric-flow residuals (must be ~0) ===')
m = sphere_mesh(12, 24)
r = solve(m, [], [1, 0, 0])
f = pressure_forces(m, r['cp'], [1, 0, 0], sref=np.pi)
print(f'sphere: CY={f["CS"]:.2e} CZ_side={f["CL"]:.2e} CD={f["CDp"]:.2e}')

w = build_wing('0012', span=6.0, chord=1.0, n_chord=16, n_span=6)
wk = build_wake_from_meta(w)
A, Brow = assemble(w, wk, kutta_mode='doublet')
r = solve(w, wk, [1, 0, 0], A=A, Brow=Brow, kutta_mode='doublet')
f = pressure_forces(w, r['cp'], [1, 0, 0], sref=6.0)
print(f'wing a=0: CL={f["CL"]:.2e} CDp={f["CDp"]:.2e} CS={f["CS"]:.2e}')

b = box(size=(1, 1, 1), n=(8, 8, 8))
r = solve(b, [], [1, 0, 0])
f = pressure_forces(b, r['cp'], [1, 0, 0], sref=1.0)
print(f'cube: CY={f["CS"]:.2e} CZ={f["CL"]:.2e} CD={f["CDp"]:.2e}')

print('=== B. Cp error map (sphere, by latitude) ===')
ms = sphere_mesh(12, 24)
rs = solve(ms, [], [1, 0, 0])
C = ms.centroid
th = np.arccos(np.clip(C[:, 0], -1, 1)) * 180 / np.pi
exact = 1 - 2.25 * np.sin(np.radians(th)) ** 2
for lo, hi in [(0, 20), (20, 40), (40, 60), (60, 80), (80, 90)]:
    sel = (th >= lo) & (th < hi)
    e = np.abs(rs['cp'][sel] - exact[sel])
    print(f'theta {lo}-{hi}: maxerr={e.max():.3f} meanerr={e.mean():.3f}')

print('=== C. conditioning ===')
for name, mm, ww in [('sphere', m, []), ('wing16x6', w, wk),
                     ('cube', b, [])]:
    A2, _ = assemble(mm, ww)
    print(f'{name}: cond={np.linalg.cond(A2):.2e} npan={mm.npanels}')

print('=== D. TE jump vs alpha (wing 24x12) ===')
w2 = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wk2 = build_wake_from_meta(w2)
A2, Brow2 = assemble(w2, wk2, kutta_mode='doublet')
nl = w2.meta['n_loop']
ns = w2.meta['n_span']
npl = nl - 1
for adeg in [0, 2, 4, 6, 8]:
    a = np.radians(adeg)
    V = [np.cos(a), 0, np.sin(a)]
    r = solve(w2, wk2, V, A=A2, Brow=Brow2, kutta_mode='doublet')
    j = ns // 2
    dte = abs(r['cp'][j * npl] - r['cp'][j * npl + npl - 1])
    print(f'alpha {adeg}: dCpTE={dte:.4f}')

print('=== E. sliver check (min area vs max Cp) ===')
from wing3d.shapes import build_f16
f16 = build_f16(scale=0.2)
print(f'f16 min area: {f16.area.min():.2e}, npan={f16.npanels}')
r = solve(f16, [], [1, 0, 0])
print(f'f16 max|Cp|: {np.abs(r["cp"]).max():.2f} (spikes if slivers dominate)')
