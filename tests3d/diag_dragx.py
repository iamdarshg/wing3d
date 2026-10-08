import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
cp = res['cp']
C = mesh.centroid
# drag direction ~ x (alpha small)
q = 0.5
vhat = V / np.linalg.norm(V)
dF = (-cp[:, None] * mesh.normal * mesh.area[:, None])
dD = (dF @ vhat) * q
dD /= 6.0  # sref
bins = np.linspace(0, 1.0, 11)
xc = C[:, 0]
print('CDp by chord station (midspan |y|<1.5 only):')
tot = 0.0
for b in range(10):
    m = (xc >= bins[b]) & (xc < bins[b + 1]) & (np.abs(C[:, 1]) < 1.5)
    print('  x=%.1f-%.1f: CDp=%+.6f' % (bins[b], bins[b + 1], dD[m].sum()))
    tot += dD[m].sum()
print('midspan-strip total: %.5f' % tot)
