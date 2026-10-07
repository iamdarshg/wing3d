import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake, assemble,
                            wing_strips_from_structured)
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12,
                  te_gap=0.001)
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
strips, _ = wing_strips_from_structured(nl, ns)
for tag, get in [('upperTE', lambda j: np.array(mesh.meta['te_points'][j]) + np.array([0, 0, 0.0005])),
                 ('mid', lambda j: np.array(mesh.meta['te_points'][j])),
                 ('lowerTE', lambda j: np.array(mesh.meta['te_points'][j]) - np.array([0, 0, 0.0005]))]:
    info = []
    for (up, lo, j) in strips:
        pj = get(j)
        pj1 = get(j + 1)
        info.append((up, lo, pj, pj1))
    wakes = build_wake(mesh, info)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print(f'{tag}: CL={f["CL"]:+.4f} CDp={f["CDp"]:+.5f}')
