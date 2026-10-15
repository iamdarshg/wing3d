"""Isolate loop-vs-profile contribution to the RAE CLI gap (M0.85)."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

from scipy.signal import savgol_filter
from tsd_m_sup import run_sup
from tsd_coupled import build_case
from tsd_campaign import build_profile, build_section
from wing3d.tsdcli import extract_section
from wing3d.geometry import read_stl, section_loop
from wing3d.tsd import section_profile_xy

mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
sec = extract_section(mesh)
yu_s = savgol_filter(sec['yu'], 41, 3)
yl_s = savgol_filter(sec['yl'], 41, 3)
prof_e = section_profile_xy(sec['x'], yu_s, yl_s)
loop_e = section_loop(sec['x'], yu_s, sec['x'], yl_s, 24)
prof_t = build_profile('rae2822')
loop_t = build_section('rae2822')

d = np.hypot(loop_e[:, 0] - loop_t[:, 0], loop_e[:, 1] - loop_t[:, 1])
print('loop max |diff| = %.5f (nose pt d=%.5f)' % (d.max(), d[24]))
print('loop_e LE pt (%.4f,%+.4f)  loop_t LE pt (%.4f,%+.4f)'
      % (loop_e[24, 0], loop_e[24, 1], loop_t[24, 0], loop_t[24, 1]))

print('== A: table loop + extracted profile')
sA = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6, section=loop_t)
LA, _, _, _ = run_sup(0.85, span=1.6, epochs=2, free_k=False,
                      profile=prof_e, section=loop_t, base=sA)
print('A CL=%.4f' % LA['CL'], flush=True)

print('== B: extracted loop + table profile')
sB = build_case(0.3, 1.25, False, True, None, 80, 32, 32, 1.6, section=loop_e)
LB, _, _, _ = run_sup(0.85, span=1.6, epochs=2, free_k=False,
                      profile=prof_t, section=loop_e, base=sB)
print('B CL=%.4f' % LB['CL'], flush=True)
