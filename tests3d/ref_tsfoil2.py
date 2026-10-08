"""Dump pyTSFoil reference: transition, Cf, edge velocity per Re."""
import os
os.environ['PYTHONUTF8'] = '1'
import numpy as np
from pytsfoil import run_airfoil_analysis

beta = np.linspace(0, np.pi, 120)
x = (1 - np.cos(beta)) / 2
yt = 5 * 0.12 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
loop = np.vstack([[[1.0, 0.0]],
                  np.stack([x[::-1], yt[::-1]], axis=1)[1:],
                  np.stack([x, -yt], axis=1)[1:]])
out = {}
for Re in [1e5, 1e6, 3e6, 1e7]:
    try:
        r = run_airfoil_analysis(loop, Mach=0.5, AoA_degrees=4.0, Re=Re,
                                 flag_IBL=True)
        ib = r.get('ibl_upper', {})
        ib2 = r.get('ibl_lower', {})
        print('Re=%.0e CL=%.4f CDf=%.5f CDtot=%.5f' % (
            Re, r.get('cl'), r.get('cd_friction'), r.get('cd_total')))
        print('  ibl_upper keys:', list(ib.keys())[:12])
        print('  ibl_lower keys:', list(ib2.keys())[:12])
        out[Re] = r
    except Exception as e:
        print('Re=%.0e FAILED: %s' % (Re, str(e)[:150]))
np.save('cases/tsfoil_ref.npy', out, allow_pickle=True)
print('saved')
