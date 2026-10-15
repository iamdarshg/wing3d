"""Cambered-imprint probe: panel + TSD CL for naca0012 vs rae2822.

Measures whether the camber-carrying panel imprint produces the missing
camber loading in the TSD field (span=1.6, alpha=1.25).
"""
import sys
import time
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from tsd_coupled import build_case
from tsd_campaign import build_section, build_profile
from tsd_m_sup import run_sup

SPAN = 1.6
for geo in ('naca0012', 'rae2822'):
    sec = build_section(geo)
    prof = build_profile(geo)
    s, res, mesh = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              SPAN, section=sec)
    muw = np.asarray(res['muw'])
    y = np.linspace(-SPAN / 2, SPAN / 2, len(muw))
    CLp = 2.0 * np.trapezoid(muw, y) / SPAN
    print('PANEL %-8s CL=%.4f (muw integral)' % (geo, CLp), flush=True)
    for M in (0.3, 0.85):
        t0 = time.time()
        L, xc, ix, j0 = run_sup(M, span=SPAN, epochs=4, free_k=False,
                                profile=prof, section=sec)
        print('TSD   %-8s M=%.2f CL=%.4f (%.0fs)' % (
            geo, M, L['CL'], time.time() - t0), flush=True)
