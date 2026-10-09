import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.shock import shock_fitted_cp, pocket_and_shock

mesh = build_wing('0012', span=1.6, chord=1.0, n_chord=24, n_span=8)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(1.25)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
C = mesh.centroid
for M in [0.85, 0.885, 0.95]:
    cs = {}
    for sgn, tag in [(1, 'up'), (-1, 'lo')]:
        m = (np.abs(C[:, 1]) < 0.2) & ((C[:, 2] * sgn) > 0.001)
        o = np.argsort(C[m, 0])
        x, cp0 = C[m, 0][o], res['cp'][m][o]
        xu = np.unique(x.round(4))
        cu = np.array([cp0[np.abs(x - xq) < 1e-3].mean() for xq in xu])
        r = shock_fitted_cp(xu, cu, M)
        pk = pocket_and_shock(xu, cu, M)
        cs[tag] = (r['x'], r['cp_corr'], r['x_shock'], pk['M1'])
    print('M=%.3f up-shock=%s M1=%.2f | lo-shock=%s M1=%.2f' % (
        M, cs['up'][2], cs['up'][3], cs['lo'][2], cs['lo'][3]))
