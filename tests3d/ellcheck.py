import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.primitives import ellipsoid
from wing3d.solver import solve

for nl, nn in [(10, 16), (16, 32), (20, 40)]:
    m = ellipsoid(1.0, 0.5, 0.5, nlat=nl, nlon=nn)
    res = solve(m, [], [1, 0, 0])
    print(nl, 'x', nn, 'npan=', m.npanels, 'maxCp=', round(float(res['cp'].max()), 3))
