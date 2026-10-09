"""Performance measurement: wing3d paths (C vs numpy) + OpenFOAM refs."""
import numpy as np
import os
import sys
import time
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

CASES = [
    ('wing-AR6 (1248 pan)', dict(code='0012', span=6.0, chord=1.0,
                                 n_chord=24, n_span=12), 6.0),
    ('wing-AR16 (576 pan)', dict(code='0012', span=1.6, chord=1.0,
                                 n_chord=24, n_span=8), 1.6),
    ('F-16-like (1600 pan)', None, 0.806),
]
print('path | case | assemble | solve | total')
for name, kw, sref in CASES:
    if kw is None:
        from wing3d.shapes import build_f16_waked
        mesh, wakes, sref, _ = build_f16_waked(scale=0.2)
    else:
        mesh = build_wing(**kw)
        wakes = build_wake_from_meta(mesh)
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    t0 = time.time()
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    tA = time.time() - t0
    t0 = time.time()
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    tS = time.time() - t0
    f = pressure_forces(mesh, res['cp'], V, sref=sref)
    tag = 'C' if not os.environ.get('WING3D_NO_C') else 'NP'
    print('%s | %s | %.1fs | %.1fs | %.1fs CL=%.4f' % (
        tag, name, tA, tS, tA + tS, f['CL']))
print('OpenFOAM refs (this host): car RANS-SST ~4min; central Euler '
      '12k ~10-19min; wing RANS 172k mesh-limited ~30min+')
