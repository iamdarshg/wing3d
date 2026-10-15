"""Panel-CL sanity test after the TE-gap fix: dfvlr (finite TE -> now
closed to 0.002) and rae (closed TE, unchanged). Panel solves only."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from tsd_coupled import build_case
from tsd_campaign import build_section

SPAN = 1.6
for geo in ('dfvlr_r4', 'rae2822', 'naca0012'):
    sec = build_section(geo)
    print('%s section TE: upper=%+.5f lower=%+.5f gap=%.5f' % (
        geo, sec[0, 1], sec[-1, 1], sec[0, 1] - sec[-1, 1]))
    s, res, mesh = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              SPAN, section=sec)
    muw = np.asarray(res['muw'])
    y = np.linspace(-SPAN / 2, SPAN / 2, len(muw))
    CLp = 2.0 * np.trapezoid(muw, y) / SPAN
    print('PANEL %-9s CL=%.4f' % (geo, CLp), flush=True)
