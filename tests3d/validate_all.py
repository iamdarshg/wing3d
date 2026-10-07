import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces, trefftz_cd
from wing3d.shapes import build_car, build_f16

print('=== NACA0012 AR6 polar (doublet-Kutta, inviscid) ===')
print('ref: LL CL(4deg)=0.33; VLM=0.39; Abbott 2D cl(4deg)=0.44')
mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
for adeg in [0, 2, 4, 6, 8]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print('alpha %d: CL=%+.4f CDp=%+.5f' % (adeg, f['CL'], f['CDp']))

print('=== Ahmed-like car, alpha=0 (inviscid, qualitative) ===')
car = build_car()
res = solve(car, [], [1, 0, 0])
f = pressure_forces(car, res['cp'], [1, 0, 0], sref=0.389 * 0.288)
print('car CD = %.4f (experiment ~0.23-0.28 with separation; '
      'inviscid omits it)' % f['CDp'])

print('=== F-16-like, alpha=0/4 (inviscid, qualitative) ===')
f16 = build_f16(scale=0.2)
for adeg in [0, 4]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(f16, [], V)
    f = pressure_forces(f16, res['cp'], V, sref=1.0)
    print('alpha %d: CL=%+.4f CDp=%+.5f' % (adeg, f['CL'], f['CDp']))
