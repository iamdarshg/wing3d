"""Wing polar validation vs public data.

References (no internet; standard published values):
- Abbott & von Doenhoff, Theory of Wing Sections: NACA0012 2D cl~0.11/deg,
  cl(4deg)~0.43-0.44, cd~0.006-0.008 at Re 3-6M.
- Prandtl lifting line, AR6 rectangular: CL(4deg) ~ 0.33.
- VLM (tests3d/vlm.py, same repo): CL(4deg) ~ 0.39.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces, trefftz_cd

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
wakes = build_wake_from_meta(mesh)
nl = mesh.meta['n_loop']
ns = mesh.meta['n_span']
npl = nl - 1
dy = 6.0 / (2 * 12)


def strip_cls(cp):
    ys, cls, ch = [], [], []
    for j in range(ns):
        ids = [j * npl + i for i in range(npl)]
        Fz = (-cp[ids, None] * mesh.normal[ids]
              * mesh.area[ids, None]).sum(axis=0)[2]
        ys.append(float(np.mean(mesh.centroid[ids][:, 1])))
        cls.append(float(Fz / dy))  # cl = L'/(q*c), q=0.5 absorbed: Fz_raw/dy
        ch.append(1.0)
    return np.array(ys), np.array(cls), np.array(ch)


print('alpha | mode      | CL     | CDp     | CDi(tr) | e     | dCpTE')
for adeg in [0, 2, 4, 6, 8]:
    a = np.radians(adeg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    for km in ['doublet', 'potential']:
        A, Brow = assemble(mesh, wakes, kutta_mode=km)
        res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode=km)
        f = pressure_forces(mesh, res['cp'], V, sref=6.0)
        ys, cls, ch = strip_cls(res['cp'])
        tr = trefftz_cd(cls, ys, ch)
        j = ns // 2
        dte = abs(res['cp'][j * npl] - res['cp'][j * npl + npl - 1])
        print(f'{adeg:5d} | {km:9s} | {f["CL"]:+.4f} | {f["CDp"]:+.5f} | '
              f'{tr["CDi"]:+.5f} | {tr["e"]:.3f} | {dte:.4f}')
