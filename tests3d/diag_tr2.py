import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.coupled import trace_streamlines
from wing3d.ibl import solve_streamline, michel_transition

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
lines, _ = trace_streamlines(mesh, res['vel'])
print('Re midspan-track transition sweep (upper tracks only):')
for ln in lines:
    if ln.get('fragment'):
        continue
    C = mesh.centroid[ln['panels']]
    if C[:, 2].mean() <= 0:
        continue
    if abs(C[:, 1].mean()) > 1.0:
        continue
    n = len(ln['s'])
    row = []
    for Re in [1e5, 3e5, 1e6, 3e6, 1e7]:
        r = solve_streamline(ln['s'], ln['Ue'], 1.0 / Re)
        row.append('%d/%d%s' % (r['i_tr'], n,
                                'B' if r.get('bubble') else ''))
    print('ymid=%+.2f Ue0=%.2f ' % (C[:, 1].mean(), ln['Ue'][0])
          + ' '.join(row))
