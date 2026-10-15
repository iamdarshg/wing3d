"""Panel baselines: CL(muw span-integral) at span=1.6 and 6.0, alpha 1.25."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')
from tsd_coupled import build_case

for span in (1.6, 6.0):
    s, res, mesh = build_case(0.3, 1.25, False, True, None, 80, 32, 32, span)
    muw = np.asarray(res['muw'])
    print('span %.1f: res keys=%s' % (span, sorted(res.keys())))
    print('  muw shape=%s sum=%.5f mean=%.5f' % (muw.shape, muw.sum(),
                                                 muw.mean()))
    # CL from the spanwise circulation: CL = 2*int(Gamma dy)/(U*S), U=1
    y = np.linspace(-span / 2, span / 2, len(muw))
    gamma = muw
    CL = 2.0 * np.trapezoid(gamma, y) / span
    print('  panel CL (span-integral of muw) = %.4f' % CL)
    # also from cp via up/lo pairing on (x,y) grid columns
    C = np.asarray(mesh.centroid)
    cp = np.asarray(res['cp'])
    upm = C[:, 2] > 0.001
    lom = C[:, 2] < -0.001
    xcoord = np.round(C[:, 0], 6)
    ycoord = np.round(C[:, 1], 6)
    pairs = {}
    for i in np.where(upm)[0]:
        pairs.setdefault((xcoord[i], ycoord[i]), [None, None])[0] = cp[i]
    for i in np.where(lom)[0]:
        pairs.setdefault((xcoord[i], ycoord[i]), [None, None])[1] = cp[i]
    tot = 0.0
    for (xx, yy), (cu, cl) in pairs.items():
        if cu is not None and cl is not None:
            tot += (cl - cu)
    # dA approx = dx*dy from counts
    dxs = np.unique(np.diff(np.unique(xcoord)))
    dys = np.unique(np.diff(np.unique(ycoord)))
    dA = (dxs[dxs > 0].min() if len(dxs) else 0.05) * (
        dys[dys > 0].min() if len(dys) else 0.1)
    print('  panel CL (cp pairing) = %.4f  (n=%d, dA~%.5f)' % (
        tot * dA / span, sum(1 for v in pairs.values()
                             if v[0] is not None and v[1] is not None), dA))
