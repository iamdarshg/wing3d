"""Transonic Mach sweep M0.8->1.1 for A330 + F-16 (and F-5E sane).

Method (simple sweep theory):
  M_n = M cos(Lambda_eff); section stays subsonic-normal up to
  M = 1/cos(Lambda) ~ 1.15 for these 30-32 deg swept wings, so KT and
  shock-fitted cp apply in the NORMAL plane.
  Section cl(M) via KT (M_n < Mcrit) or shock_fitted_cp (M_n >= Mcrit).
  Config CL(M) = CL(M0.3panel) * [cl_sec(M)/cl_sec(M0.3)] preserving
  config AR/fuse character while carrying section transonic physics.
  M_n >= 1 (not reached) would need Ackeret; flagged if it does.
Error bands: KT +-7%, shockfit +-15%, Ackeret +-20% (flagged).
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shapes import build_f16_waked, build_a330, build_f5
from wing3d.solver import solve, assemble
from wing3d.forces import pressure_forces
from wing3d.transonic import karman_tsien_cp, m_crit
from wing3d.shock import shock_fitted_cp

CFGS = {
    'A330': (build_a330, dict(scale=0.02), 30.0, 0.33),
    'F-16': (build_f16_waked, dict(scale=0.2), 32.0, 0.19),
    'F-5E': (build_f5, dict(scale=0.1), 24.0, 0.33),
}
MS = [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1]
ADE = 2.0  # sweep at 2 deg (transonic cruise-like)


def grab_section(mesh, cp, y0, dy=0.05):
    """Exposed mid-semispan section -> (x, cpu, cpl)."""
    C = mesh.centroid
    N = mesh.normal
    out = {}
    for sgn, tag in [(1, 'up'), (-1, 'lo')]:
        m = (np.abs(C[:, 1] - y0) < dy) & ((N[:, 2] * sgn) > 0.3)
        if m.sum() < 6:
            return None
        o = np.argsort(C[m, 0])
        x, c = C[m, 0][o], np.asarray(cp)[m][o]
        xn = (x - x.min()) / max(x.max() - x.min(), 1e-9)
        xu = np.unique(xn.round(3))
        cu = np.array([c[np.abs(xn - xq) < 5e-3].mean() for xq in xu])
        out[tag] = (xu, cu)
    return out


def section_cl(sec, M, Mn, mc):
    """Section cl at normal Mach with method flag."""
    xg = np.linspace(0, 1, 400)
    if Mn < 1.0:
        if Mn < mc:
            u = np.interp(xg, sec['up'][0], sec['up'][1])
            l = np.interp(xg, sec['lo'][0], sec['lo'][1])
            u = karman_tsien_cp(u, Mn)
            l = karman_tsien_cp(l, Mn)
            method = 'KT'
        else:
            ru = shock_fitted_cp(sec['up'][0], sec['up'][1], Mn)
            rl = shock_fitted_cp(sec['lo'][0], sec['lo'][1], Mn)
            u = np.interp(xg, ru['x'], ru['cp_corr'])
            l = np.interp(xg, rl['x'], rl['cp_corr'])
            method = 'shockfit'
        cl = float(np.trapezoid(l - u, xg))
        return cl, method
    # Mn >= 1 (supersonic normal): no linear model (transonic similarity
    # regime); needs TSD/OF. Flag instead of Ackeret (which blows up at
    # beta->0).  Higher-Mn valid range reached only by low-sweep wings.
    a = np.radians(ADE)
    return float('nan'), 'Euler/TSD-needed'


for name, (b, kw, lam_deg, y0) in CFGS.items():
    mesh, wakes, sref, info = b(**kw)
    C = mesh.centroid
    span = C[:, 1].max() - C[:, 1].min()
    AR = span ** 2 / sref
    skip = mesh.meta.get('skipped', [])
    a = np.radians(ADE)
    V = np.array([np.cos(a), 0, np.sin(a)])
    A_, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    res = solve(mesh, wakes, V, A=A_, Brow=Brow, kutta_mode='doublet')
    f0 = pressure_forces(mesh, res['cp'], V, sref=sref, skip=skip)
    CL0 = f0['CL']  # M0 baseline (incompressible panel)
    sec = grab_section(mesh, res['cp'], y0, 0.05)
    if sec is None:
        print('=== %s: section not found at y=%.2f, abort ===' % (
            name, y0), flush=True)
        continue
    mc = m_crit(min(sec['up'][1].min(), sec['lo'][1].min()))
    cl03, _ = section_cl(sec, 0.3, 0.3, mc)  # baseline section
    cosl = np.cos(np.radians(lam_deg))
    inv_c = 1.0 / cosl  # M where M_n = 1
    print('=== %s: AR=%.2f sweep=%.0fdeg section-y=%.2f Mcrit_norm=%.2f '
          'M_sonic_norm=%.2f CL0=%.3f cl0section=%.3f ===' % (
              name, AR, lam_deg, y0, mc, inv_c, CL0, cl03), flush=True)
    print('  M  |  Mn  | sec-cl | conf-CL | shock_x | method')
    for M in MS:
        Mn = M * cosl
        cl, method = section_cl(sec, M, Mn, mc)
        confCL = CL0 * cl / cl03 if np.isfinite(cl) else float('nan')
        shx = None
        if method == 'shockfit':
            ru = shock_fitted_cp(sec['up'][0], sec['up'][1], Mn)
            shx = ru['x_shock']
        band = {'KT': 7, 'shockfit': 15}.get(method, 25)
        flag = '' if Mn < 1.0 else ' (Mn>=1: no linear model)'
        print('  %.2f | %.2f | %+.3f | %+.3f | %s | %s (+-%d%%)%s' % (
            M, Mn, cl, confCL,
            ('%.2f' % shx) if shx else '-',
            method, band, flag), flush=True)