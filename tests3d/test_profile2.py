import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
import numpy as np
from tsd_campaign import build_profile
from gen_transonic_stl import rae2822_xy

# raw coords
xu, yu, xl, yl = rae2822_xy()
print('raw: xu[0..3]=%s yu[0..3]=%s' % (np.round(xu[:3], 3), np.round(yu[:3], 4)))
print('raw: xl[0..3]=%s yl[0..3]=%s' % (np.round(xl[:3], 3), np.round(yl[:3], 4)))

p = build_profile('rae2822')
xg = np.linspace(0, 1, 400)
hu = p['half'](xg)
mean = p['mean'](xg)
for xq in [0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98]:
    i = int(np.argmin(np.abs(xg - xq)))
    print('x=%.2f half=%+.4f mean=%+.5f' % (xg[i], hu[i], mean[i]))