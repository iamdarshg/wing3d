"""Mach + Re sweeps v2: section-based transonic, Cf-based viscous errors.

Part A (Mach, a=2/4): config CL@M0.3 vs LL+PG (error%); wing-section
  cl vs M via KT below Mcrit / shock-fitted above (uncertainty bands).
Part B (Re, a=2): coupled; Cf_wet = CDf*S/Awet vs Prandtl-Schlichting
  (true viscous error%); CDtot reported with config-fidelity caveat.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f5, build_a330
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces
from wing3d.transonic import karman_tsien_cp, m_crit
from wing3d.shock import shock_fitted_cp
from wing3d.coupled import solve_coupled

CFGS = {
    'F-5E': (build_f5, dict(scale=0.1), 0.30),
    'A330': (build_a330, dict(scale=0.02), 0.33),
}
MS = [0.3, 0.5, 0.65, 0.8, 0.85, 0.9, 0.95]
RES = [1e5, 3e5, 1e6, 3e6, 1e7]


def section_curves(mesh, cp, y0, dy=0.06):
    C = mesh.centroid
    N = mesh.normal
    out = {}
    for sgn, tag in [(1, 'up'), (-1, 'lo')]:
        m = (np.abs(C[:, 1] - y0) < dy) & ((N[:, 2] * sgn) > 0.3)
        o = np.argsort(C[m, 0])
        x, c = C[m, 0][o], np.asarray(cp)[m][o]
        xn = (x - x.min()) / max(x.max() - x.min(), 1e-9)
        xu = np.unique(xn.round(3))
        cu = np.array([c[np.abs(xn - xq) < 5e-3].mean() for xq in xu])
        out[tag] = (xu, cu)
    return out


for name, (builder, kw, y0) in CFGS.items():
    mesh, wakes, sref, info = builder(**kw)
    C = mesh.centroid
    span = C[:, 1].max() - C[:, 1].min()
    AR = span ** 2 / sref
    Awet = float(mesh.area.sum())
    print('=== %s: npan=%d S=%.4f AR=%.2f Awet/S=%.2f ===' % (
        name, mesh.npanels, sref, AR, Awet / sref), flush=True)
    for adeg in [2.0, 4.0]:
        a = np.radians(adeg)
        V = np.array([np.cos(a), 0, np.sin(a)])
        A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
        res = solve(mesh, wakes, V, A=A, Brow=Brow,
                    kutta_mode='doublet')
        f0 = pressure_forces(mesh, res['cp'], V, sref=sref)
        beta = np.sqrt(1 - 0.09)
        ll = 2 * np.pi * AR / (AR + 2) * a / beta
        print('%s a=%.0f: CL(M0.3)=%.3f vs LLPG=%.3f err=%+.1f%%' % (
            name, adeg, f0['CL'],
            ll, (f0['CL'] - ll) / ll * 100), flush=True)
        sec = section_curves(mesh, res['cp'], y0)
        mc = m_crit(min(sec['up'][1].min(), sec['lo'][1].min()))
        print('  section Mcrit=%.2f' % mc, flush=True)
        print('  M | sec-cl | shock_x | method')
        for M in MS:
            if M < mc:
                xg = np.linspace(0, 1, 200)
                cu = karman_tsien_cp(sec['up'][1], M)
                cl_ = karman_tsien_cp(sec['lo'][1], M)
                cl = float(np.trapezoid(
                    np.interp(xg, sec['lo'][0], cl_) -
                    np.interp(xg, sec['up'][0], cu), xg))
                print('  %.2f | %.3f | - | KT(+-7%%)' % (M, cl),
                      flush=True)
            elif M < 1.0:
                r = shock_fitted_cp(sec['up'][0], sec['up'][1], M)
                rl = shock_fitted_cp(sec['lo'][0], sec['lo'][1], M)
                xg = np.linspace(0, 1, 200)
                cl = float(np.trapezoid(
                    np.interp(xg, rl['x'], rl['cp_corr']) -
                    np.interp(xg, r['x'], r['cp_corr']), xg))
                print('  %.2f | %.3f | %s | shockfit(+-15%%)' % (
                    M, cl, r['x_shock']), flush=True)
    print('-- Re sweep (a=2): Cf error vs Prandtl-Schlichting --')
    print('Re | CL | CDtot | CDf | Cf_wet | Cf_PS | err% | it')
    a = np.radians(2.0)
    V = np.array([np.cos(a), 0, np.sin(a)])
    Lref = sref / span
    for Re in RES:
        res, f, inf = solve_coupled(mesh, wakes, V, Re=Re, Lref=Lref,
                                    itmax=8, relax=0.2)
        h = inf['history'][-1]
        CL = h['CL'] / sref
        CD = h['CD'] / sref
        CDf = CD - f['CDp'] / sref
        Cf = CDf * sref / Awet
        ReL = Re  # Re based on Lref=mean chord already
        CfPS = 0.455 / (np.log10(ReL) ** 2.58)
        print('%.0e | %.4f | %.5f | %.5f | %.5f | %.5f | %+.1f%% | %d'
              % (Re, CL, CD, CDf, Cf, CfPS,
                 (Cf - CfPS) / CfPS * 100, inf['iterations']),
              flush=True)
