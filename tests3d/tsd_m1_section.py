"""E4: same-config RCA -- TSD M1.0 CENTERLINE section vs OF slab anchor.

OF slab (walls, span_in 1.2) ~= 2D section; TSD AR6 full-wing CL carries
tip relief. Centerline section cl is the apples-to-apples number.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_m1 import run_m1

s = run_m1(alpha_deg=1.25)
L = s.fp_loads()
jy = np.asarray(s.wing_jy)
j0 = jy[np.argmin(np.abs(s.yc[jy]))]
jtip = jy[-1]
dxi = s.dx
ix = np.asarray(s.wing_ix)
xc = s.xc[ix] / s.chord


def section(j):
    d = (L['cpl'][ix, j] - L['cpu'][ix, j]) * dxi[:len(ix)]
    return float(d.sum() / s.chord)


print('full-wing CL = %.4f (OF slab anchor 0.0787)' % L['CL'])
print('centerline section cl = %.4f' % section(j0))
print('tip section cl = %.4f' % section(jtip))
# shock location at centerline: max adverse jump in cpu
cpu = L['cpu'][ix, j0]
dj = np.diff(cpu)
ish = int(np.argmax(dj))
print('centerline shock x~=%.2f (OF 0.80)' % xc[ish])
print('centerline Cpmin=%.2f' % L['cpu'][ix, j0].min())
np.save('cases/openfoam/transonic/tsd_m1_sec.npy',
        np.stack([xc, L['cpu'][ix, j0], L['cpl'][ix, j0]], axis=1))
print('saved section')
