import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import build_wake_from_meta
from wing3d.coupled import solve_coupled

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=16, n_span=6)
wakes = build_wake_from_meta(mesh)
a = np.radians(4)
res, f, info = solve_coupled(mesh, wakes, [np.cos(a), 0, np.sin(a)],
                             Re=3e6, Lref=1.0, itmax=20, relax=0.2)
print('iters:', info['iterations'],
      'sep frac:', round(info['bl']['sep_fraction'], 3))
for h in info['history']:
    print({k: round(v, 5) if isinstance(v, float) else v
           for k, v in h.items()})
print('TRUE CL ~ history.CL/6, TRUE CD ~ history.CD/6')
