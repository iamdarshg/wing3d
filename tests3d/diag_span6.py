"""Extra diagnostic: TSD naca0012 at span=6.0 (near-2D) vs span=1.6 at M0.85.

Tests the 'OF slab is wall-bound / quasi-2D' hypothesis for the symmetric
case: if OF naca (0.0728) is near the span-6 TSD value, the slab's y-walls
(which the wing SPANS exactly) suppress tip relief in the symmetric case.
"""
import sys
import time
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from tsd_coupled import build_case
from tsd_campaign import build_section, build_profile
from tsd_m_sup import run_sup

for span in (6.0,):
    sec = build_section('naca0012')
    prof = build_profile('naca0012')
    s, res, mesh = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              span, section=sec)
    muw = np.asarray(res['muw'])
    y = np.linspace(-span / 2, span / 2, len(muw))
    CLp = 2.0 * np.trapezoid(muw, y) / span
    print('PANEL naca span=%.1f CL=%.4f' % (span, CLp), flush=True)
    t0 = time.time()
    L, xc, ix, j0 = run_sup(0.85, span=span, epochs=4, free_k=False,
                            profile=prof, section=sec)
    print('TSD   naca M0.85 span=%.1f CL=%.4f (%.0fs)' % (
        span, L['CL'], time.time() - t0), flush=True)
