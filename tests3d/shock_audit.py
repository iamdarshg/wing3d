"""Shock-region audit across the Mach sweep.

For each aircraft section and each M (0.7..1.1):
  - supersonic-pocket extent (start/end, PG vs Cp* sonic-line crossing)
  - pocket islands / mixed-stream diagnostics (num pockets, gaps,
    subsonic-supersonic interleaving near the sonic line)
  - shock location + jump + M1, wave drag Cdw and its share of total Cd
  - near-boundary checks: pocket touches LE (x=0) or TE (x=1), shock
    at/after TE, wake-side behavior
  - ERROR % of shock_x and CDw vs the OpenFOAM anchors where they exist
    (straight NACA0012 incompressible section -> shock_fitted vs OF).

Anchor shock_x the user wants to compare against (alpha=1.25deg):
    M0.8 0.35 | M0.885 0.40 | M0.95 0.60 | M1.0 0.80 | M1.1 0.70
OF section CD (total, incl ~0.04 coarse floor): use CDw_ref = CD-CD0.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shock import (shock_fitted_cp, shock_region_metrics,
                          pocket_and_shock, mach_from_cp)
from wing3d.transonic import karman_tsien_cp, m_crit
from wing3d.shapes import build_f16_waked, build_a330

MACHES = [0.7, 0.75, 0.8, 0.85, 0.885, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15]
# OF anchors (M, alpha=1.25): (CL, CD, shock_x, cpmin_up)
OF = {(0.8, 1.25): (0.0686, 0.0485, 0.35, -1.27),
      (0.885, 1.25): (0.0696, 0.0512, 0.40, -0.74),
      (0.95, 1.25): (0.0830, 0.0650, 0.60, -0.76),
      (1.0, 1.25): (0.0790, 0.0770, 0.80, -0.72),
      (1.1, 1.25): (0.0533, 0.1030, 0.70, -0.60)}
OF_CD0 = 0.044  # coarse-mesh floor proxy at M0.8 a1.25


def load_section(tag):
    d = np.load('tmp_sec0012_%s.npy' % tag)
    return d[0], d[1]


def pocket_islands(pk, x):
    """Split the sup bool mask into separate supersonic islands."""
    sup = pk['pocket']
    if not np.any(sup):
        return []
    edges = np.diff(sup.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0]
    if len(starts) == 0 or len(ends) == 0:
        return [(float(x[0]), float(x[-1]))]
    if starts[0] > ends[0]:
        starts = np.r_[0, starts]
    if ends[-1] < starts[-1]:
        ends = np.r_[ends, len(sup) - 1]
    return [(float(x[s]), float(x[min(e, len(x) - 1)])) for s, e in zip(starts, ends)]


def audit_section(x0, cp0, name, a_wake_cut=1.0):
    xg = np.linspace(0, 1, 300)
    u0, l0 = load_section('up'), load_section('lo')
    xu, cu = u0[0], u0[1]
    xl, cl = l0[0], l0[1]
    print('=' * 78)
    print('%s (straight NACA0012 alpha=1.25deg)' % name)
    print('%-5s %-5s %-7s %-7s %-8s %-7s %-7s %-7s %-8s %-8s %s' % (
        'M', 'Mcri', 'pocket', 'pkt_%', 'mix_%', 'shock', 'M1',
        'Cd_wav', 'shock%', 'Cdw%', 'islands/bndry'))
    tot_cd_waves = {}
    for M in MACHES:
        mc = m_crit(cu.min())
        if M >= 1.0:
            # PG singular above sonic -> report surrogate/OF only
            of = OF.get((round(M, 3), 1.25)) or OF.get((round(M, 3) - 0.05, 1.25))
            if of is None:
                continue
            print('%-5.2f %-5s %-7s %-7s %-8s %-7s %-7s %-7s %-8s %-8s %s' % (
                M, '--', '--', '--', '--', '%.2f' % of[2], '--',
                '--', '--', '--',
                'no-lin-model (needs TSD/OF)'))
            continue
        # KT below Mcrit (no pocket), shock-fitted above
        if M < mc:
            cu_ = karman_tsien_cp(np.interp(xg, xu, cu), M)
            cl_ = karman_tsien_cp(np.interp(xg, xl, cl), M)
            pk = pocket_and_shock(xg, cu_, M, smooth=0)
            has_pocket = pk['x_shock'] is not None
        else:
            r = shock_fitted_cp(xg, np.interp(xg, xu, cu), M)
            has_pocket = r['x_shock'] is not None
        if not has_pocket:
            print('%-5.2f %-5.2f %-7s %-7s %-8s %-7.3s %-7s %-7s %-8s %-8s %s' % (
                M, mc, 'none', '0', '0', '-', '-', '-', '-', '-',
                'subcritical (KT)'))
            continue
        # re-run with the ORIGINAL incompressible cp (not KT-amplified)
        # because pocket logic wants Cp0 / beta internally
        r = shock_fitted_cp(xg, np.interp(xg, xu, cu), M)
        met = shock_region_metrics(xg, np.interp(xg, xu, cu), M)
        islands = pocket_islands(met, xg)
        bnd = ''
        if met['x_start'] is not None and met['x_start'] < 0.01:
            bnd += 'LE-touch '
        if met['x_end'] is not None and met['x_end'] > 0.99:
            bnd += 'TE-touch '
        if r['x_shock'] is not None and r['x_shock'] >= 1.0:
            bnd += 'shock>=TE '
        if len(islands) > 1:
            bnd += 'islands:%d ' % len(islands)
        # wave drag share of total Cd (estimate: Cdw + 0.008 profile floor)
        cd_tot_est = met['cdw'] + 0.008
        # vs OF anchors
        key = (round(M, 3), 1.25)
        of = OF.get(key)
        serr = werr = None
        if of is not None and M < 1.0:
            of_shock = of[2]
            if of_shock:
                serr = (r['x_shock'] - of_shock) / of_shock * 100
            of_cdw = max(of[1] - OF_CD0, 0.0001)
            werr = (met['cdw'] - of_cdw) / of_cdw * 100
        print('%-5.2f %-5.2f %-7s %-7.3f %-8.2f %-7.2f %-7.2f %-8.4f %-8.1f %-8.1f %s' % (
            M, mc,
            '%.2f-%.2f' % ((met['x_start'] or 0), (met['x_end'] or 1)),
            met['pocket_fraction'], met['mixed_fraction'] * 100,
            r['x_shock'] if r['x_shock'] is not None else float('nan'),
            met['M1'], met['cdw'],
            serr if serr is not None else float('nan'),
            werr if werr is not None else float('nan'),
            bnd if bnd else 'ok'))
        tot_cd_waves[M] = met['cdw']
    return tot_cd_waves


x0 = np.linspace(0, 1, 24)
audit_section(x0, None, 'REFERENCE section')
print()
print('(shock%%/Cdw%% = error vs OpenFOAM anchor at that M;')
print(' Cdw_ref = OF CD - %.3f coarse floor)' % OF_CD0)