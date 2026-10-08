import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.coupled import trace_streamlines, smooth_surface_field
from wing3d.ibl import solve_streamline

d = np.load('cases/tsfoil_ref.npy', allow_pickle=True).item()
u = d[3e6]['ibl_upper']
ts = u['s']
# tsfoil s is arc length from LE? check range
print('tsfoil s range:', ts.min().round(3), ts.max().round(3),
      'x_tr:', u['x_tr'], 'Cfmax:', u['cf'].max().round(4))
# their cf vs x: need x mapping; xx_foil available?
r = d[3e6]
print('keys xx:', 'xx' in r, 'xx_foil' in r)
if 'xx_foil' in r:
    xx = np.asarray(r['xx_foil'])
    print('xx_foil[:5]:', xx[:5].round(3) if xx.ndim == 1 else xx.shape)

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
vm = np.linalg.norm(res['vel'], axis=1)
vs = smooth_surface_field(mesh, vm, passes=1)
dirn = res['vel'] / np.maximum(vm, 1e-12)[:, None]
lines, _ = trace_streamlines(mesh, dirn * vs[:, None])
best = None
for ln in lines:
    if ln.get('fragment'):
        continue
    C = mesh.centroid[ln['panels']]
    if C[:, 2].mean() > 0.01 and abs(C[:, 1].mean()) < 0.5:
        best = ln
        break
rr = solve_streamline(best['s'], best['Ue'], 1.0 / 3e6)
print('ours i_tr=%d/%d bubble=%s Cfmax=%.4f Cfavg=%.5f' % (
    rr['i_tr'], len(best['s']), rr.get('bubble'), rr['cf'].max(),
    rr['cf'].mean()))
# Cf profile comparison at matched stations
n = len(best['s'])
ii = np.linspace(0, n - 1, 9).astype(int)
print('ours s :', np.round(best['s'][ii], 3))
print('ours Cf:', np.round(rr['cf'][ii], 5))
cf_ts = u['cf']
jj = np.linspace(0, len(cf_ts) - 1, 9).astype(int)
print('tsfoil Cf:', np.round(cf_ts[jj], 5))
