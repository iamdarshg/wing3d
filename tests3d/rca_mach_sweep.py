"""RCA for the F-5E/F-16/A330 Mach-sweep 23% error + F-16 instability.

Findings (2026-10-11):
1. KT + Cp<=-4 clip + skipped-mask: F-5E 0.155 vs LLPG 0.151 = +2.9% (GOOD).
   A330 0.182 vs 0.186 = -2.1% (GOOD).
2. F-16: no TE wake + tip slivers (taper on clipped delta) =>
   KT amplification of unclipped root/tip suction => CL 0.167 vs LLPG
   0.153 (+8.7%) but Cd=0.0011 (spurious pressure-drag cancellation) =>
   not sweepable.
3. The 23% number was measured with a WRONG pipeline (raw cp, no KT,
   no mask). Real pipeline errors are <=9% at M0.3.
Fixes implemented:
- F-16 rebuild with wing TE wake + tip taper smoothing.
- forces skip-mask always applied (buried panels excluded from KT and
  pressure integration).
- KT bounded by section Cp* (no unbounded amplification at peak suction).
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f5, build_f16_waked, build_a330
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces
from wing3d.transonic import karman_tsien_cp

a = np.radians(2.0)
V = np.array([np.cos(a), 0, np.sin(a)])
for name, b, kw in [('F-5E', build_f5, dict(scale=0.1)),
                    ('F-16', build_f16_waked, dict(scale=0.2)),
                    ('A330', build_a330, dict(scale=0.02))]:
    mesh, wakes, sref, info = b(**kw)
    skip = mesh.meta.get('skipped', [])
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    cp = karman_tsien_cp(np.maximum(res['cp'], -4.0), 0.3)
    f = pressure_forces(mesh, cp, V, sref=sref, skip=skip)
    span = info['span']
    AR = span ** 2 / sref
    ll = 2 * np.pi * AR / (AR + 2) / np.sqrt(1 - 0.09) * np.sin(a)
    print('%-5s CL=%.3f LLPG=%.3f err=%+.1f%% Cd=%.4f' % (
        name, f['CL'], ll, (f['CL'] - ll) / ll * 100, f['CDp']))
