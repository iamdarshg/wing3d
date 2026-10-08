import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.coupled import trace_streamlines
from wing3d.ibl import solve_streamline

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
lines, _ = trace_streamlines(mesh, res['vel'])
ups = [ln for ln in lines if not ln.get('fragment')
       and mesh.centroid[ln['panels']][:, 2].mean() > 0.01
       and abs(mesh.centroid[ln['panels']][:, 1].mean()) < 1.0]
print('upper mid tracks:', len(ups))
for ln in ups[:4]:
    row = []
    for Re in [1e5, 1e6, 3e6, 1e7]:
        r = solve_streamline(ln['s'], ln['Ue'], 1.0 / Re)
        # transition x/c approx: s fraction * chord-ish
        row.append('%.2f%s' % (r['i_tr'] / len(ln['s']),
                               'B' if r.get('bubble') else ''))
    print('ymid=%+.2f ' % (mesh.centroid[ln['panels']][:, 1].mean())
          + ' '.join(row) + ' Cf3e6=%.5f' % solve_streamline(
              ln['s'], ln['Ue'], 1.0 / 3e6)['cf'].mean())
