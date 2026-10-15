"""Post-SavGol slope comparison: extraction vs table profile (RAE)."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

from scipy.signal import savgol_filter
from tsd_campaign import build_profile
from wing3d.tsdcli import extract_section
from wing3d.geometry import read_stl
from wing3d.tsd import section_profile_xy

mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
sec = extract_section(mesh)
prof_t = build_profile('rae2822')

variants = {'raw': (sec['yu'], sec['yl'])}
for w in (21, 41, 61):
    variants['sg%d' % w] = (savgol_filter(sec['yu'], w, 3),
                            savgol_filter(sec['yl'], w, 3))

xg = np.linspace(0, 1, 4001)
print('%-8s %9s | %s' % ('x', 'tbl sq', ' '.join(
    '%-17s' % v for v in variants)))
for xc in [0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.76, 0.85, 0.9, 0.95, 0.99]:
    row = ['%-8.4f %+9.4f ' % (xc, prof_t['camber'](xc))]
    for v, (yu, yl) in variants.items():
        p = section_profile_xy(sec['x'], yu, yl)
        row.append('| %+7.4f %+7.4f ' % (p['camber'](xc), p['slope'](xc)))
    print(''.join(row))

print()
print('%-8s %9s | %s' % ('x', 'tbl st', ' '.join(
    '%-17s' % v for v in variants)))
for xc in [0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.76, 0.85, 0.9, 0.95, 0.99]:
    row = ['%-8.4f %+9.4f ' % (xc, prof_t['slope'](xc))]
    for v, (yu, yl) in variants.items():
        p = section_profile_xy(sec['x'], yu, yl)
        row.append('| %+7.4f %+7.4f ' % (p['slope'](xc), p['half'](xc)))
    print(''.join(row))

print()
for v, (yu, yl) in variants.items():
    p = section_profile_xy(sec['x'], yu, yl)
    xg2 = np.linspace(0, 1, 2001)
    cb = p['camber'](xg2)
    print('%s: max slope_c=%.4f   camber p2p via %.4f' %
          (v, cb.max(), 0))
