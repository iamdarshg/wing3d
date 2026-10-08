import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
C = mesh.centroid
amax = mesh.area.max()
print('area min/max: %.2e %.4f' % (mesh.area.min(), amax))
for thr in [0.0, 1e-7, 1e-6, 1e-5, 1e-4]:
    skip = np.where(mesh.area < thr)[0] if thr > 0 else None
    f = pressure_forces(mesh, res['cp'], V, sref=6.0, skip=skip)
    nskip = 0 if skip is None else len(skip)
    print('area_thr=%.0e skip=%d CL=%.4f CDp=%.5f maxCp_kept=%.1f' % (
        thr, nskip, f['CL'], f['CDp'],
        np.abs(res['cp'] if skip is None else np.delete(res['cp'], skip)).max()))
