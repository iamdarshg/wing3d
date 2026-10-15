"""Print the exact _wing_bc inputs (profile slope/camber) for both paths."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

from tsd_campaign import build_profile
from wing3d.tsdcli import extract_section
from wing3d.geometry import read_stl
from wing3d.tsd import section_profile_xy

mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
sec_e = extract_section(mesh)
prof_e = section_profile_xy(sec_e['x'], sec_e['yu'], sec_e['yl'])
prof_t = build_profile('rae2822')

print('nose values: yu[0]=%+.5f yl[0]=%+.5f  (extraction)'
      % (sec_e['yu'][0], sec_e['yl'][0]))
xg = np.linspace(0, 1, 2000)
sl_t = prof_t['slope'](xg)
cb_t = prof_t['camber'](xg)
sl_e = prof_e['slope'](xg)
cb_e = prof_e['camber'](xg)
fc_t = 0.5 * (prof_t['slope'](xg) + 0j).real  # placeholder no-op
print('%-8s %10s %10s | %10s %10s | %9s %9s' %
      ('x', 'slope_t T', 'slope_c T', 'slope_t E', 'slope_c E',
       'd_slope_t', 'd_slope_c'))
for xc in [0.0, 0.00125, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03,
           0.05, 0.1, 0.2, 0.5, 0.7, 0.76, 0.9, 0.99]:
    i = np.argmin(np.abs(xg - xc))
    print('%-8.5f %+10.4f %+10.4f | %+10.4f %+10.4f | %+9.4f %+9.4f'
          % (xg[i], sl_t[i], cb_t[i], sl_e[i], cb_e[i],
             sl_e[i] - sl_t[i], cb_e[i] - cb_t[i]))
print('max |d_slope_t| = %.4f  max |d_slope_c| = %.4f'
      % (np.abs(sl_e - sl_t).max(), np.abs(cb_e - cb_t).max()))
