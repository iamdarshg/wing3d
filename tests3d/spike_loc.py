import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble

for gap in [0.0, 0.002]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12,
                      te_gap=gap)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    i = int(np.argmin(res['cp']))
    print('gap', gap, 'Cpmin=%+.2f' % res['cp'].min(), 'at x=%.3f z=%+.4f'
          % (mesh.centroid[i][0], mesh.centroid[i][2]))
