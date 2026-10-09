import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f5
from wing3d.solver import solve, assemble
from wing3d.shock import shock_fitted_cp

mesh, wakes, sref, info = build_f5(scale=0.1)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
C = mesh.centroid
N = mesh.normal
# mid-semispan wing section (normal-based sides for the low wing)
for adeg in [2.0, 4.0]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    for M in [0.9, 0.95, 1.0, 1.1, 1.25]:
        cs = {}
        for sgn, tag in [(1, 'up'), (-1, 'lo')]:
            m = (np.abs(C[:, 1] - 0.30) < 0.06) & ((N[:, 2] * sgn) > 0.3)
            if m.sum() < 8:
                cs[tag] = (np.array([0.0, 1.0]), np.array([0.0, 0.0]),
                           None)
                continue
            o = np.argsort(C[m, 0])
            x, cp0 = C[m, 0][o], res['cp'][m][o]
            xn = (x - x.min()) / max(x.max() - x.min(), 1e-9)
            xu = np.unique(xn.round(3))
            cu = np.array([cp0[np.abs(xn - xq) < 5e-3].mean()
                           for xq in xu])
            if M < 1.0:
                r = shock_fitted_cp(xu, cu, M)
                cs[tag] = (r['x'], r['cp_corr'], r['x_shock'])
            else:
                # supersonic: Ackeret section estimate instead
                beta = np.sqrt(max(M * M - 1, 1e-9))
                cla = 4 * a / beta
                cs[tag] = (None, None, None)
        if M < 1.0:
            xg = np.linspace(0, 1, 200)
            cpu = np.interp(xg, cs['up'][0], cs['up'][1])
            cpl = np.interp(xg, cs['lo'][0], cs['lo'][1])
            cl = float(np.trapezoid(cpl - cpu, xg))
            print('a=%.1f M=%.2f section-cl=%.4f sh_up=%s' % (
                adeg, M, cl, cs['up'][2]))
        elif M == 1.0:
            print('a=%.1f M=1.00 no linear estimate (transonic '
                  'similarity regime; needs OpenFOAM/TSD)' % adeg)
        else:
            beta = np.sqrt(max(M * M - 1, 1e-9))
            print('a=%.1f M=%.2f Ackeret-2D cl=%.4f (thin-wing est)' % (
                adeg, M, 4 * a / beta))
