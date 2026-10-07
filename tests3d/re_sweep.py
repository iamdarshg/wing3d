"""Error vs Re (speed sweep): viscous-coupled wing across Re decades.

Truth anchors (public data):
- Abbott/von Doenhoff NACA0012 2D drag bucket + Ladson data.
- OpenFOAM SST (car validated; wing section via 2D when available).
Finds the most problematic Re for the IBL/turbulence modeling.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import build_wake_from_meta
from wing3d.coupled import solve_coupled

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
a = np.radians(4)
V = [np.cos(a), 0, np.sin(a)]
print('Re | CL | CDtot | CDf | sep_frac | trans_x/c (mid) | iters')
for Re in [1e5, 3e5, 1e6, 3e6, 1e7]:
    res, f, info = solve_coupled(mesh, wakes, V, Re=Re, Lref=1.0,
                                 itmax=12, relax=0.2)
    bl = info['bl']
    # true coefficients (history sref=1 -> divide by 6)
    h = info['history'][-1]
    print(f'{Re:.0e} | {h["CL"] / 6:.4f} | {h["CD"] / 6:.5f} | '
          f'sep={bl["sep_fraction"]:.3f} sep_lines={bl["sep_lines"]}/'
          f'{bl["n_lines"]} it={info["iterations"]}', flush=True)
