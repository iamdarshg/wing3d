import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.shock import cp_star

# truth: (span, nc, ns, alpha, M, xshock_up, xshock_lo)
CASES = [
    (12.0, 32, 8, 2.0, 0.75, 0.52, None),
    (1.6, 24, 8, 1.25, 0.8, 0.35, None),
    (1.6, 24, 8, 1.25, 0.885, 0.40, 0.40),
    (1.6, 24, 8, 1.25, 0.95, 0.60, 0.60),
]
cfgs = {}
for span, nc, ns, adeg, M, xu, xl in CASES:
    key = (span, nc, ns, adeg)
    if key not in cfgs:
        mesh = build_wing('0012', span=span, chord=1.0, n_chord=nc,
                          n_span=ns, cosine=(span > 2))
        wakes = build_wake_from_meta(mesh)
        A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
        a = np.radians(adeg)
        V = np.array([np.cos(a), 0, np.sin(a)])
        res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
        cfgs[key] = (mesh, res, V)
    mesh, res, V = cfgs[key]
    C = mesh.centroid
    m = (np.abs(C[:, 1]) < 0.5) & (C[:, 2] > 0.001)
    o = np.argsort(C[m, 0])
    x, cp0 = C[m, 0][o], res['cp'][m][o]
    X = np.unique(x.round(4))
    CU = np.array([cp0[np.abs(x - xq) < 1e-3].mean() for xq in X])
    b = max(np.sqrt(1 - M * M), 1e-9)
    cppg = CU / b
    cps = cp_star(M)
    sup = np.where(cppg < cps)[0]
    ipk = sup[np.argmin(cppg[sup])]
    cmin = cppg[ipk]
    print('M=%.3f pocket %.3f-%.3f' % (M, X[sup[0]], X[sup[-1]]))
    for fr in [0.3, 0.5, 0.7]:
        tgt = cmin + fr * (cps - cmin)
        # first x past peak where Cp recovers above tgt
        cands = np.where((np.arange(len(X)) > ipk) & (cppg > tgt))[0]
        xs = X[cands[0]] if len(cands) else -1
        print('   recov %.1f -> xsh=%.3f (truth %s)' % (fr, xs, xu))
