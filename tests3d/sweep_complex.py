"""Mach + Re sweeps: F-5E and A330 (A320-class airliner surrogate).

Part A (Mach, alpha=0/2/4): one panel solve per alpha, KT correction per M.
  Error = |CL_KT - CL_LLPG| / CL_LLPG (lifting-line + PG theory).
Part B (Re, alpha=2): full viscous-coupled per Re.
  CD0_est = CD - CL^2/(pi*AR*e); error vs public CD0 figures.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f5, build_a330
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces
from wing3d.transonic import karman_tsien_cp, kt_status
from wing3d.coupled import solve_coupled

CFGS = {
    'F-5E': (build_f5, dict(scale=0.1), 0.0215, 'Jane s CD0~0.020-0.023'),
    'A330': (build_a330, dict(scale=0.02), 0.0200, 'cruise CD0~0.019-0.021'),
}
MS = [0.3, 0.5, 0.65, 0.8, 0.85, 0.9]
RES = [1e5, 3e5, 1e6, 3e6, 1e7]

for name, (builder, kw, cd0_pub, note) in CFGS.items():
    mesh, wakes, sref, info = builder(**kw)
    C = mesh.centroid
    span = C[:, 1].max() - C[:, 1].min()
    AR = span ** 2 / sref
    print('=== %s: npan=%d S=%.4f span=%.3f AR=%.2f ===' % (
        name, mesh.npanels, sref, span, AR), flush=True)
    # ---- Part A: Mach ----
    print('-- Mach sweep --')
    print(' a\\M | ' + ' | '.join('%.2f' % m for m in MS) +
          ' || LLPG-err%')
    for adeg in [0.0, 2.0, 4.0]:
        a = np.radians(adeg)
        V = np.array([np.cos(a), 0, np.sin(a)])
        A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
        res = solve(mesh, wakes, V, A=A, Brow=Brow,
                    kutta_mode='doublet')
        cp0 = np.maximum(res['cp'], -4.0)  # clip LE singularity (proven path)
        row = []
        errs = []
        for M in MS:
            cp = karman_tsien_cp(cp0, M)
            f = pressure_forces(mesh, cp, V, sref=sref)
            st = kt_status(M, cp0.min())
            cla_ll = 2 * np.pi * AR / (AR + 2) * a / np.sqrt(
                max(1 - M ** 2, 1e-9))
            err = abs(f['CL'] - cla_ll) / max(abs(cla_ll), 1e-9) * 100
            errs.append(err)
            flag = '' if st == 'ok' else '*'
            row.append('%.3f%s' % (f['CL'], flag))
        print(' %3.0f  | %s || %.1f%%' % (
            adeg, ' | '.join(row), np.mean(errs)), flush=True)
    # ---- Part B: Re ----
    print('-- Re sweep (a=2) --')
    print('Re | CL | CDtot | CDp | CDf | CD0est | CD0err% | it')
    a = np.radians(2.0)
    V = np.array([np.cos(a), 0, np.sin(a)])
    Lref = sref / span  # mean chord
    for Re in RES:
        res, f, inf = solve_coupled(mesh, wakes, V, Re=Re, Lref=Lref,
                                    itmax=8, relax=0.2)
        h = inf['history'][-1]
        CL = h['CL'] / sref
        CD = h['CD'] / sref
        CD0 = CD - CL ** 2 / (np.pi * AR * 0.85)
        err = (CD0 - cd0_pub) / cd0_pub * 100
        print('%.0e | %.4f | %.5f | %.5f | %.5f | %.5f | %+.1f%% | %d' % (
            Re, CL, CD, f['CDp'] / sref,
            (CD - f['CDp'] / sref), CD0, err, inf['iterations']),
            flush=True)
    print('public: %s' % note, flush=True)
print('(* = KT flagged; LLPG = lifting-line+PG theory)')
