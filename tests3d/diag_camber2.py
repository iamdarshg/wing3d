"""Cambered-imprint probe #2: panel + TSD CL for dfvlr_r4 (expect ~OF 0.267
at M0.85) plus biconvex sanity at M0.85."""
import sys
import time
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from tsd_coupled import build_case
from tsd_campaign import build_section, build_profile
from tsd_m_sup import run_sup

SPAN = 1.6
for geo, machs in (('dfvlr_r4', (0.30, 0.85)), ('biconvex', (0.85,))):
    sec = build_section(geo)
    prof = build_profile(geo)
    s, res, mesh = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              SPAN, section=sec)
    muw = np.asarray(res['muw'])
    y = np.linspace(-SPAN / 2, SPAN / 2, len(muw))
    CLp = 2.0 * np.trapezoid(muw, y) / SPAN
    print('PANEL %-9s CL=%.4f (muw integral)' % (geo, CLp), flush=True)
    for M in machs:
        t0 = time.time()
        L, xc, ix, j0 = run_sup(M, span=SPAN, epochs=4, free_k=False,
                                profile=prof, section=sec)
        print('TSD   %-9s M=%.2f CL=%.4f (%.0fs)' % (
            geo, M, L['CL'], time.time() - t0), flush=True)
