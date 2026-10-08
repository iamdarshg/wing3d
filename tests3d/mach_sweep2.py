"""Speed (Mach) sweep: panel + Karman-Tsien vs Prandtl-Glauert theory."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble,
                           karman_tsien)
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
f0 = pressure_forces(mesh, res['cp'], V, sref=6.0)
print('M    CL KT      CL PG-theory  err')
for M in [0.0, 0.2, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75]:
    cp = karman_tsien(res['cp'], M) if M > 0 else res['cp']
    f = pressure_forces(mesh, cp, V, sref=6.0)
    pg = f0['CL'] / max(np.sqrt(1 - M ** 2), 1e-9)
    err = (f['CL'] - pg) / pg * 100 if M > 0 else 0.0
    print('%.2f %.4f     %.4f        %+.1f%%' % (M, f['CL'], pg, err),
          flush=True)
