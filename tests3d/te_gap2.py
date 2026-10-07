import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for gap in [0.0, 0.001, 0.002, 0.005]:
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12,
                      te_gap=gap)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    print(f'te_gap={gap}: CL={f["CL"]:+.4f} CDp={f["CDp"]:+.5f} '
          f'npan={mesh.npanels}')
print('(ref LL 0.33, VLM 0.39)')
