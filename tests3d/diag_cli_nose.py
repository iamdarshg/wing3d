"""Quick check of the cos-clustered extraction nose."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from scipy.signal import savgol_filter
from wing3d.tsdcli import extract_section
from wing3d.geometry import read_stl, section_loop
from tsd_campaign import build_section

mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
sec = extract_section(mesh)
print('x [:6] :', np.array2string(sec['x'][:6], precision=5))
print('yu[:6] :', np.array2string(sec['yu'][:6], precision=5))
print('yl[:6] :', np.array2string(sec['yl'][:6], precision=5))
print('x [-4:]:', np.array2string(sec['x'][-4:], precision=5))
print('yu[-4:]:', np.array2string(sec['yu'][-4:], precision=5))
print('yl[-4:]:', np.array2string(sec['yl'][-4:], precision=5))
print('camber=%.4f thick=%.4f' % (sec['camber_max'], sec['thick_max']))

yu_s = savgol_filter(sec['yu'], 41, 3)
yl_s = savgol_filter(sec['yl'], 41, 3)
yu_s[0], yu_s[-1] = sec['yu'][0], sec['yu'][-1]
yl_s[0], yl_s[-1] = sec['yl'][0], sec['yl'][-1]
loop_e = section_loop(sec['x'], yu_s, sec['x'], yl_s, 24)
loop_t = build_section('rae2822')
d = np.hypot(loop_e[:, 0] - loop_t[:, 0], loop_e[:, 1] - loop_t[:, 1])
print('loop max |diff|=%.5f  (nose d=%.5f, TE d=%.5f, TE2 d=%.5f)'
      % (d.max(), d[24], d[0], d[48]))
print('loop_e LE (%.4f,%+.4f)  loop_t LE (%.4f,%+.4f)'
      % (loop_e[24, 0], loop_e[24, 1], loop_t[24, 0], loop_t[24, 1]))
