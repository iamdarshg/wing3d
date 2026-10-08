"""Complex-geometry campaign: F-16 / A-330 / Shuttle AoA sweeps (timed).

Inviscid panel + wakes; Trefftz induced drag as internal truth.
Published numbers are approximate (config differs: no strakes vortex,
no separation, no friction) -- compared with explicit caveats.
"""
import numpy as np
import sys
import time
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces, trefftz_cd
from wing3d.shapes import build_f16_waked, build_a330, build_shuttle
from wing3d.forces import spanwise_loading


def run_cfg(name, builder, kw, alphas, sref_note=''):
    t0 = time.time()
    mesh, wakes, sref, info = builder(**kw)
    t_build = time.time() - t0
    # verify merge preserved panel order/counts
    t0 = time.time()
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    t_asm = time.time() - t0
    print('== %s: panels=%d wakes=%d sref=%.3f build=%.1fs asm=%.1fs' %
          (name, mesh.npanels, len(wakes), sref, t_build, t_asm))
    out = []
    skip = mesh.meta.get('skipped', None)
    nskip = 0 if skip is None else len(np.atleast_1d(skip))
    print('  masked buried panels: %d' % nskip)
    for adeg in alphas:
        a = np.radians(adeg)
        V = np.array([np.cos(a), 0, np.sin(a)])
        t0 = time.time()
        res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
        t_sol = time.time() - t0
        f = pressure_forces(mesh, res['cp'], V, sref=sref, skip=skip)
        out.append((adeg, f['CL'], f['CDp'], t_sol))
        print('  a=%4.1f CL=%+.4f CDp=%.5f sol=%.1fs' % (adeg, f['CL'],
                                                         f['CDp'], t_sol))
    return out


print('--- F-16-like (subsonic, V~200m/s class) ---')
print('ref(approx, clean subsonic): CD0~0.02-0.03, CLa~0.08/deg')
f16 = run_cfg('f16', build_f16_waked, dict(scale=0.2), [0, 4, 8, 12])
print('--- A-330-like (cruise CL~0.5, L/D~19-20 published) ---')
a330 = run_cfg('a330', build_a330, dict(scale=0.05), [0, 2, 4])
print('--- Shuttle-orbiter-like (approach, pub L/Dmax~4-5) ---')
shu = run_cfg('shuttle', build_shuttle, dict(scale=0.05), [0, 5, 10, 15])
