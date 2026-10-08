"""Reference runs: pyTSFoil (NASA TSFOIL2 + IBL) on NACA0012.

Same architecture as wing3d (TSD + Thwaites/Michel/Head IBL) -> direct
comparison of transition, Cf, drag.
"""
import os
os.environ['PYTHONUTF8'] = '1'
import numpy as np
from pytsfoil import run_airfoil_analysis

# NACA0012 loop (x, y), TE -> upper -> LE -> lower -> TE
beta = np.linspace(0, np.pi, 120)
x = (1 - np.cos(beta)) / 2
yt = 5 * 0.12 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
xu = x[::-1]
yu = yt[::-1]
xl = x
yl = -yt
loop = np.vstack([[[1.0, 0.0]],
                  np.stack([xu, yu], axis=1)[1:],
                  np.stack([xl, yl], axis=1)[1:]])  # TE->up->LE->lo->TE
print('npts:', len(loop))
for Re in [1e5, 1e6, 3e6, 1e7]:
    try:
        out = run_airfoil_analysis(loop, Mach=0.5, AoA_degrees=4.0, Re=Re,
                                   flag_IBL=True)
        keys = list(out.keys())
        print('Re=%.0e keys:' % Re, keys)
        for k in ['CL', 'CD', 'CDf', 'CDp', 'xtr_u', 'xtr_l', 'transition',
                  'x_tr_upper', 'x_tr_lower']:
            if k in out:
                print('  ', k, '=', out[k])
    except Exception as e:
        print('Re=%.0e FAILED: %s' % (Re, str(e)[:200]))
