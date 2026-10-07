import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake, assemble,
                            wing_strips_from_structured)
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
meta = mesh.meta
nl, ns = meta['n_loop'], meta['n_span']
strips, _ = wing_strips_from_structured(nl, ns)
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
for gap in [0.0, 0.002, 0.005, 0.01]:
    info = []
    for (up, lo, j) in strips:
        pa = np.array(meta['te_points'][j]) + np.array([gap, 0, 0])
        pb = np.array(meta['te_points'][j + 1]) + np.array([gap, 0, 0])
        info.append((up, lo, pa, pb))
    wakes = build_wake(mesh, info)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print(f'wake inset {gap}: CL={f["CL"]:+.4f} CDp={f["CDp"]:+.5f} '
          f'Cpmin={res["cp"].min():+.2f}')
