"""Isolate the farfield effect on the free-circulation solves.

All previous A-ladder rows used WING3D_FARFIELD=wake (the default in
kutta mode), whose box values come from the *panel* wake muw that this
repo already documents as ill-conditioned (CLp +0.63/-7.50/-1.28 at
span 1.6/3.2/6.0, M0.85 a1.25). If that box is what kills the lift, the
kutta rows are not the culprit and the farfield has to be rebuilt from
Gamma. This script runs {farfield off, wake} x {mode 0, 1, 3} at M0.3,
2 epochs, symmetric section, and prints CL side by side.

Usage: python tsd_a9_ff.py
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

CASES = [('0', 'off'), ('0', 'wake'), ('3', 'off'), ('3', 'wake')]


def run(mode, ff, M=0.3, epochs=2):
    os.environ['WING3D_KUTTA'] = mode
    os.environ['WING3D_FARFIELD'] = ff
    os.environ['WING3D_IMPRINT'] = 'll'
    os.environ['WING3D_ENTROPY'] = '0'
    os.environ['WING3D_BL'] = '0'
    from tsd_campaign import build_section
    from tsd_coupled import build_case
    from tsd_m_sup import run_sup
    sec = build_section('dfvlr_r4')
    base, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6,
                              section=sec)
    L, xc, ix, j0 = run_sup(M, span=1.6, epochs=epochs, free_k=False,
                            profile=None, section=sec, base=base)
    rec = {'mode': mode, 'farfield': ff, 'CL': float(L['CL'])}
    if 'kutta' in L:
        rec['Gamma_mid'] = float(L['kutta']['Gamma'][len(L['kutta']['Gamma'])
                                                   // 2])
        rec['hist'] = [(h['ep'], round(h['CL'], 4), round(h['gamma_mid'], 4))
                       for h in L['kutta']['hist']]
    print('FARFIELD %s mode%s CL=%+.5f %s' % (
        ff, mode, rec['CL'], rec.get('hist', '')), flush=True)
    return rec


if __name__ == '__main__':
    out = [run(m, f) for m, f in CASES]
    p = os.path.join(REPO, 'tests3d', 'a9_farfield_probe.json')
    json.dump(out, open(p, 'w'), indent=1)
    print('saved', p)
