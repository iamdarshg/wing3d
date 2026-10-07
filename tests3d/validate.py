import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for nc, ns in [(16, 6), (24, 10), (32, 14)]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=nc, n_span=ns)
    wakes = build_wake_from_meta(mesh)
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    for km in ['doublet', 'potential']:
        A, Brow = assemble(mesh, wakes, kutta_mode=km)
        res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode=km)
        f = pressure_forces(mesh, res['cp'], V, sref=6.0)
        print(f'{nc}x{ns} {km:9s} CL={f["CL"]:+.4f} CDp={f["CDp"]:+.5f}')
print('(lifting-line 0.33, VLM 0.39)')
