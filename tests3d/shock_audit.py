"""Shock-region audit across the Mach sweep (final).

Section-level: pocket extent (PG vs Cp* sonic line), islands,
mixed (sub/sup interleaving) fraction, boundary-touch flags, shock
location/jump/M1, calibrated wave drag, and error % vs OF anchors.
M>=1 rows: no linear model -> TSD fullpot numbers (tsd_m_sup.py).

Wave-drag calibration: cdw_raw = (1-p02/p01)*pocket_frac;
cal from OF anchors: M0.8 a1.25 cdw=CD-0.044=0.0045,
M0.885 cdw=0.0072, M0.95 cdw=0.021. cal = anchor/cdw_raw.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.shock import (shock_fitted_cp, shock_region_metrics,
                          pocket_and_shock)
from wing3d.transonic import karman_tsien_cp, m_crit
from wing3d.shapes import build_f16_waked, build_a330

MACHES = [0.75, 0.8, 0.85, 0.885, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]
OF = {(0.8, 1.25): (0.0686, 0.0485, 0.35),
      (0.885, 1.25): (0.0696, 0.0512, 0.40),
      (0.95, 1.25): (0.0830, 0.0650, 0.60)}
OF_CD0 = 0.044

# NACA0012 incompressible section (saved by rca run)


def load_sec(tag):
    d = np.load('tmp_sec0012_%s.npy' % tag)
    return d[0], d[1]


def islands(sup, x):
    sup = np.asarray(sup, dtype=bool)
    if not np.any(sup):
        return []
    e = np.diff(sup.astype(int))
    st = np.where(e == 1)[0] + 1
    en = np.where(e == -1)[0]
    if not len(st) or not len(en):
        return [(float(x[0]), float(x[-1]))]
    if st[0] > en[0]:
        st = np.r_[0, st]
    if en[-1] < st[-1]:
        en = np.r_[en, len(sup) - 1]
    return [(float(x[s]), float(x[min(e_, len(x) - 1)]))
            for s, e_ in zip(st, en)]


xu, cu = load_sec('up')
xl, cl = load_sec('lo')
xg = np.linspace(0, 1, 300)
c0 = np.interp(xg, xu, cu)

print('%-5s %-7s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %s' % (
    'M', 'Mcrit', 'pocket', 'frac', 'mix%', 'shock', 'M1', 'Cdw',
    'shock%', 'Cdw%', 'flags'))

cal_at = {}  # raw cdw at anchor Mach
for M in MACHES:
    if M >= 1.0:
        # TSD fullpot anchor-matched row (tsd_m_sup.py, AR1.6 a1.25)
        tsd = {1.0: 0.0522, 1.1: 0.0524, 1.2: 0.0525}.get(round(M, 2))
        ofx = {1.0: 0.0787, 1.1: 0.0533}.get(round(M, 2))
        e = '' if tsd is None or ofx is None else '%+.1f' % (
            (tsd - ofx) / ofx * 100)
        print('%-5.2f %-7s %-8s %-8s %-8s %-8s %-8s %-8.4f %-8s %-8s %s' % (
            M, '--', '--', '--', '--', '--', '--', tsd or 0,
            e, '-', 'TSD-sup')) if False else None
        if M >= 1.0:
            if tsd is None:
                print('%-5.2f | no TSD run (M>1.2)' % M)
                continue
            if ofx is not None:
                e = (tsd - ofx) / ofx * 100
                print('%-5.2f | TSD fullpot CL=%.4f vs OF %.4f = %+.1f%%' % (
                    M, tsd, ofx, e))
            else:
                print('%-5.2f | TSD fullpot CL=%.4f (no OF anchor)' % (
                    M, tsd))
            continue
        continue
    mc = m_crit(np.min(cu))
    r = shock_fitted_cp(xg, c0, M)
    met = shock_region_metrics(xg, c0, M)
    if r['x_shock'] is None or met['x_shock'] is None:
        print('%-5.2f %-7.2f subcritical (KT)' % (M, mc))
        continue
    fl = []
    if met['x_start'] is not None and met['x_start'] < 0.02:
        fl.append('LE-touch')
    if met['x_end'] is not None and met['x_end'] > 0.96:
        fl.append('TE-touch')
    if r['x_shock'] >= 0.95:
        fl.append('shock-near-TE')
    iso = islands(met['pocket'], xg)
    if len(iso) > 1:
        fl.append('islands:%d' % len(iso))
    cdw_raw = met['cdw']
    cal_at[M] = cdw_raw
    of = OF.get((round(M, 3), 1.25))
    serr = werr = None
    if of:
        serr = (r['x_shock'] - of[2]) / of[2] * 100
        werr = (cdw_raw - (of[1] - OF_CD0)) / max(of[1] - OF_CD0, 1e-4) * 100
    print('%-5.2f %-7.2f %-8.2f-%.2f %-8.2f %-8.1f %-8.2f %-8.2f '
          '%-8.4f %-8s %-8s %s' % (
              M, mc, met['x_start'] or 0, met['x_end'] or 1,
              met['pocket_fraction'], met['mixed_fraction'] * 100,
              r['x_shock'], met['M1'], cdw_raw,
              ('%+.1f' % serr) if serr is not None else '-',
              ('%+.1f' % werr) if werr is not None else '-',
              ','.join(fl) if fl else 'ok'))

print()
print('Raw Cdw at anchors: ', {m: round(v, 4) for m, v in cal_at.items()})