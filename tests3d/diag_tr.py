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
print('nlines:', len(lines))
ups = los = fr = 0
for ln in lines:
    if ln.get('fragment'):
        fr += 1
        continue
    zmean = mesh.centroid[ln['panels']][:, 2].mean()
    if zmean > 0:
        ups += 1
    else:
        los += 1
print('upper tracks:', ups, 'lower tracks:', los, 'fragments:', fr)
# transition per full track at Re 3e6 + where Michel crosses
for ln in lines:
    if ln.get('fragment'):
        continue
    zmean = mesh.centroid[ln['panels']][:, 2].mean()
    r = solve_streamline(ln['s'], ln['Ue'], 1.0 / 3e6)
    n = len(ln['s'])
    print('side=%+d len=%d i_tr=%d (%.2f) bubble=%s Cfavg=%.5f' % (
        1 if zmean > 0 else -1, n, r['i_tr'], r['i_tr'] / n,
        r.get('bubble'), r['cf'].mean()))
    break
