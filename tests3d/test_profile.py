import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
import numpy as np
from tsd_campaign import build_profile

for geo in ['naca0012', 'rae2822', 'dfvlr_r4', 'biconvex']:
    p = build_profile(geo)
    xg = np.linspace(0.02, 0.98, 30)
    half = p['half'](xg)
    cam = p['camber'](xg)
    print('%-10s t/c=%.3f camber_max=%.4f camber_range=[%+.3f, %+.3f]' % (
        geo, 2 * half.max(), np.abs(cam).max(), cam.min(), cam.max()))