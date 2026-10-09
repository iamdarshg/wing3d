import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.shock import shock_fitted_cp

# fine uniform input (no cosine slivers, better LE)
mesh = build_wing('0012', span=1.6, chord=1.0, n_chord=48, n_span=6,
                  cosine=False)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(1.25)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
C = mesh.centroid
of = {0.85: 0.0706, 0.885: 0.0696, 0.95: 0.083}
errs = []
for M, clo in sorted(of.items()):
    cs = {}
    for sgn, tag in [(1, 'up'), (-1, 'lo')]:
        m = (np.abs(C[:, 1]) < 0.2) & ((C[:, 2] * sgn) > 0.001)
        o = np.argsort(C[m, 0])
        x, cp0 = C[m, 0][o], res['cp'][m][o]
        xu = np.unique(x.round(4))
        cu = np.array([cp0[np.abs(x - xq) < 1e-3].mean() for xq in xu])
        r = shock_fitted_cp(xu, cu, M)
        cs[tag] = (r['x'], r['cp_corr'], r['x_shock'])
    xg = np.linspace(0, 1, 401)
    cpu = np.interp(xg, cs['up'][0], cs['up'][1])
    cpl = np.interp(xg, cs['lo'][0], cs['lo'][1])
    cl = float(np.trapezoid(cpl - cpu, xg))
    e = (cl - clo) / clo * 100
    errs.append(abs(e))
    print('%.3f CL=%.4f vs %.4f %+.1f%% (sh %s)' % (M, cl, clo, e,
                                                    cs['up'][2]))
print('mean abs err (fine input): %.1f%%' % (sum(errs) / len(errs)))
