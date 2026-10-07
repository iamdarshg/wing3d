import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.primitives import transform
from wing3d.solver import assemble, solve

ht = build_wing('0012', span=5.0 * 0.2, chord=1.8 * 0.2, taper=0.5,
                sweep_deg=40.0, n_chord=10, n_span=6, cosine=False)
ht = transform(ht, translate=(0.88 * 3.0, 0, 0.05 * 0.2))
A, Brow = assemble(ht, [])
print('htail alone cond:', round(float(np.linalg.cond(A)), 1))
res = solve(ht, [], [1, 0, 0])
print('htail alone max|Cp|:', round(float(np.abs(res['cp']).max()), 2))
