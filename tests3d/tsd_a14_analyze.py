"""A14 tables: read every tests3d/a14_*.json and print the three
measurements the round is about.

  python tsd_a14_analyze.py [pattern]
"""
import glob
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'


def row(p):
    d = json.load(open(p))
    g = d.get('gridinfo', {})
    r1 = np.asarray(d.get('resid_phys_mid_alpha', [np.nan]))
    r2 = np.asarray(d.get('resid_phys2_mid_alpha', [np.nan]))
    # direct normalisation measure: the field's downwash response divided by
    # the kernel's for the SAME doublet sheet, over the stations where the
    # kernel is not small (|we|>0.25 alpha). This is the "1.36x" number.
    ratio = float('nan')
    if 'w_field' in d:
        wf0 = np.asarray(d['w_field'])
        j = wf0.shape[1] // 2
        al = np.radians(d['alpha_deg'])
        wf = wf0[:, j] / al
        we = np.asarray(d['w_eq'])[:, j] / al
        m = (np.abs(we) > 0.25) & (np.abs(wf) > 1e-9)
        if m.sum():
            ratio = float(np.median(wf[m] / we[m]))
    return dict(
        f=os.path.basename(p), M=d['M'], grid=d['grid'],
        thick=d.get('thick', '0'), bl=d.get('bl', '0'),
        ep=d['epochs'], CL=d['CL'],
        err=d.get('err_pct', float('nan')),
        gm=d.get('gamma_mid', float('nan')),
        tip=d.get('gamma_tip_ratio', float('nan')),
        r1m=float(r1.mean()), r1a=float(np.abs(r1).max()),
        r2m=float(r2.mean()) if 'w_thick' in d else float('nan'),
        r2a=float(np.abs(r2).max()) if 'w_thick' in d else float('nan'),
        nf=d.get('norm_factor', float('nan')),
        zu=g.get('stencil_z_up', float('nan')),
        zl=g.get('stencil_z_lo', float('nan')),
        dzb=g.get('dz_below_k0', float('nan')),
        dza=g.get('dz_above_k0', float('nan')),
        nc=g.get('n_chord', nan0()), ns=g.get('n_span', nan0()),
        ratio=ratio,
        prof=d.get('profile'), secs=d.get('secs', 0.0))


def nan0():
    return float('nan')


def main(pat='tests3d/a14_*.json'):
    fs = sorted(glob.glob(os.path.join(REPO, pat)))
    if not fs:
        print('no files match', pat)
        return
    hdr = ('%-46s %6s %-11s %s %s %2s %8s %7s %7s %6s %6s %6s %6s' % (
        'file', 'M', 'grid', 't', 'b', 'ep', 'CL', 'err%', 'G_mid',
        'tip/mid', 'res1/a', 'res2/a', 'wf/we'))
    print(hdr)
    print('-' * len(hdr))
    for p in fs:
        r = row(p)
        print('%-46s %6.3f %-11s %s %s %2d %8.5f %+7.2f %7.4f %6.3f '
              '%+6.3f %+6.3f %6.2f' % (
                  r['f'][:46], r['M'], r['grid'], r['thick'], r['bl'],
                  r['ep'], r['CL'], r['err'], r['gm'], r['tip'],
                  r['r1m'], r['r2m'], r['ratio']))
    print('\nres2/a = mean resid_phys AFTER the reference subtraction; '
          'res1/a = before.  t = WING3D_KUTTA_THICK, b = WING3D_BL.\n')
    grid = [row(p) for p in fs if row(p)['M'] == 0.3 and row(p)['prof'] is False]
    if grid:
        print('%-11s %5s %5s %9s %9s %9s %8s %7s %7s' % (
            'grid', 'n_ch', 'n_sp', 'z_up/c', 'z_lo/c', 'dz_below',
            'res1/a', 'normfac', 'wf/we'))
        for r in grid:
            print('%-11s %5.0f %5.0f %9.5f %9.5f %9.4f %+8.3f %7.3f %7.2f' % (
                r['grid'], r['nc'], r['ns'], r['zu'], r['zl'], r['dzb'],
                r['r1m'], r['nf'], r['ratio']))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'tests3d/a14_*.json')
