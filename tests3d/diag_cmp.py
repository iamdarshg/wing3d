import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.coupled import trace_streamlines

d = np.load('cases/tsfoil_ref.npy', allow_pickle=True).item()
r = d[3e6]
u = r['ibl_upper']
print('TSFOIL Re3e6 upper: xtr=%.3f' % u['x_tr'])
s, ue, th, H, cf = u['s'], u['ue'], u['theta'], u['H'], u['cf']
# ue here: edge velocity IN WHAT UNITS? (TSD scaled?) print raw + normalized
print('tsfoil s[:6]:', np.round(s[:6], 4))
print('tsfoil ue[:6]:', np.round(ue[:6], 4))
print('tsfoil ue max:', ue.max().round(3), 'at s=', s[np.argmax(ue)].round(3))
print('tsfoil theta[0], mid, tr:', th[0], th[len(th)//2], th[u['i_tr']])
print('tsfoil H at tr:', H[u['i_tr']].round(2), ' Cf avg:', cf.mean().round(5))

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
lines, _ = trace_streamlines(mesh, res['vel'])
best = None
for ln in lines:
    if ln.get('fragment'):
        continue
    C = mesh.centroid[ln['panels']]
    if C[:, 2].mean() > 0.01 and abs(C[:, 1].mean()) < 0.5:
        best = ln
        break
print('ours Ue max:', best['Ue'].max().round(3), ' Cf-less Ue[:6]:',
      np.round(best['Ue'][:6], 3))
print('ours s[:6]:', np.round(best['s'][:6], 4))
