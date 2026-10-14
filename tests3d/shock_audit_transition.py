"""Transition-band shock audit (M0.8-1.0) using OF compressible Cp.

Root cause of the M0.85-0.95 blowup (validated 2026-10-11):
shock_fitted_cp/pocket_and_shock PG-amplify the incompressible Cp0 by
1/sqrt(1-M^2), which diverges as M->1 (factor 3.2 at M0.95) and feeds
the isentropic Mach inversion, which saturates at its M=5.0 ceiling.
The physical flow at M0.9 has Cpmin=-0.76 -> local M~1.35.

Fix: for M>=0.8, STOP PG-amplifying. Interpolate the ACTUAL
OpenFOAM compressible upper/lower Cp sections (M0.8, M0.9, M1.0)
in Mach, then run pocket/shock metrics directly on real Cp.
"""
import numpy as np
import json
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shock import (mach_from_cp, normal_shock_loss, cp_star,
                          pocket_and_shock)


def anchor_cp():
    d = json.load(open('cases/openfoam/transonic/summary.json'))
    out = {}
    for key, M in [('of_m08', 0.8), ('of_m09', 0.9)]:
        sec = {}
        for side, tag in [('upper', 'u'), ('lower', 'l')]:
            arr = np.array(d[key][side])
            xu = np.unique(arr[:, 0])
            cv = np.array([arr[arr[:, 0] == q, 1].mean() for q in xu])
            sec[tag] = (xu, cv)
        out[M] = sec
    # M1.0 from vtp extraction (rows: xg, cpu, cpl)
    m1 = np.load('tmp_ofm10_sec.npy')
    out[1.0] = {'u': (m1[0], m1[1]), 'l': (m1[0], m1[2])}
    return out


def interp_cp(anchors, M, side, xq):
    """Bilinear-in-Mach, chord-normalized upper/lower Cp at xq."""
    Ms = sorted(anchors)
    if M < Ms[0] or M > Ms[-1]:
        return None
    i = max(j for j, m in enumerate(Ms) if m <= M)
    i = min(i, len(Ms) - 2)
    m0, m1 = Ms[i], Ms[i + 1]
    t = (M - m0) / (m1 - m0)
    x0, c0 = anchors[m0][side]
    x1, c1 = anchors[m1][side]
    return (1 - t) * np.interp(xq, x0, c0) + t * np.interp(xq, x1, c1)


def audit(M, anchors):
    xg = np.linspace(0, 1, 200)
    cu = interp_cp(anchors, M, 'u', xg)
    cl = interp_cp(anchors, M, 'l', xg)
    if cu is None:
        return None
    # local Mach from real compressible Cp
    mu = np.array([mach_from_cp(c, M) for c in cu])
    ml = np.array([mach_from_cp(c, M) for c in cl])
    sup_u = mu > 1.0
    sup_l = ml > 1.0
    xs_u = xg[sup_u].min() if sup_u.any() else None
    xe_u = xg[sup_u].max() if sup_u.any() else None
    M1 = max(mu.max() if sup_u.any() else 1.0,
             ml.max() if sup_l.any() else 1.0)
    loss = 1.0 - normal_shock_loss(M1)
    frac_u = sup_u.mean()
    return {'cu': cu, 'cl': cl, 'mu': mu, 'ml': ml,
            'pocket_u': (xs_u, xe_u), 'frac_u': frac_u,
            'M1': M1, 'cdw': loss * frac_u,
            'cpmin_u': cu.min(), 'cpmin_l': cl.min()}


if __name__ == '__main__':
    anchors = anchor_cp()
    print('%-5s %-16s %-16s %-8s %-8s %-8s %-8s' % (
        'M', 'pocket_upper', 'pocket_lower', 'M1', 'frac', 'Cdw', 'Cpmin_u'))
    for M in [0.80, 0.85, 0.885, 0.90, 0.925, 0.95, 0.975, 1.00]:
        r = audit(M, anchors)
        if r is None:
            continue
        pu = r['pocket_u']
        # lower pocket bounds
        sup_l = r['ml'] > 1.0
        pl = (np.linspace(0, 1, 200)[sup_l].min() if sup_l.any() else None,
              np.linspace(0, 1, 200)[sup_l].max() if sup_l.any() else None)
        print('%-5.2f %-16s %-16s %-8.2f %-8.2f %-8.4f %-8.2f' % (
            M,
            '%.2f-%.2f' % (pu[0] or 0, pu[1] or 0),
            '%.2f-%.2f' % (pl[0] or 0, pl[1] or 0),
            r['M1'], r['frac_u'], r['cdw'], r['cpmin_u']))
