import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
import numpy as np
from tsd_m_sup import run_sup
from tsd_campaign import build_profile

for geo in ['naca0012', 'rae2822']:
    prof = build_profile(geo)
    L, xc, ix, j0 = run_sup(1.1, span=1.6, epochs=6, free_k=False,
                             profile=prof)
    print('%s M1.1: CL=%.4f (camber now modeled)' % (geo, L['CL']),
          flush=True)