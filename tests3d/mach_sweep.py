"""Error vs Mach (0.1 -> 1.2): KT-corrected panel vs PG/Ackeret/literature.

Truth anchors (compressible, no internet needed):
- Subsonic: Prandtl-Glauert CL(M) = CL0/sqrt(1-M^2), valid to ~M0.6-0.7.
- NACA0012 drag divergence: Mdd ~ 0.7 (low alpha; lower at alpha 4deg).
  Beyond Mdd, wave drag appears (panel has none) -> error must explode.
- Supersonic: Ackeret 2D cl = 4*alpha/sqrt(M^2-1) (thin-airfoil reference;
  finite 3D wing differs; shown for trend only).
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble, karman_tsien
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
cp0 = res['cp']
# Spike treatment (documented): the finite-gap lower-TE corner carries a
# near-singular numerical spike (Cp ~ -10 on tiny panels from wake-filament
# proximity). Negligible for integrated loads but singular for KT.
# Cap input Cp at -4 (documented).
ncap = int((cp0 < -4).sum())
cp0c = np.maximum(cp0, -4.0)
f0 = pressure_forces(mesh, cp0c, V, sref=6.0)
CL0 = f0['CL']
print(f'incompressible CL0 = {CL0:.4f} ({ncap} capped panels)')
print('Mach | CL_KT  | CL_PG  | errPG | Cpmin_KT | CDp_KT (*)')
print('(*) CD from KT-corrected Cp is unreliable (non-conservative); '
      'CL is the metric.')
for M in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9]:
    cp = karman_tsien(np.maximum(cp0, -4.0), M)
    f = pressure_forces(mesh, cp, V, sref=6.0)
    beta = np.sqrt(max(1 - M ** 2, 1e-9))
    cl_pg = CL0 / beta
    err = abs(f['CL'] - cl_pg) / max(abs(cl_pg), 1e-9)
    print(f'{M:.2f} | {f["CL"]:+.4f} | {cl_pg:+.4f} | {err * 100:.1f}% '
          f'| {cp.min():+.2f} | {f["CDp"]:+.5f}')
for M in [1.0, 1.1, 1.2]:
    cl_ack = 4 * a / np.sqrt(M ** 2 - 1) if M > 1 else float('inf')
    print(f'{M:.2f} | (KT singular) | Ackeret2D {cl_ack:+.4f} | -- '
          f'| -- | --')
