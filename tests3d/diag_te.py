import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from tsd_recirc import build_case
from scipy.optimize import newton_krylov

s = build_case(0.3, 4.0, True)
nphi = s.nx * s.ny * s.nz


def func(x):
    s.phi[:] = x.reshape(s.nx, s.ny, s.nz)
    R = s.residual()
    R[0, :, :] = 0
    R[-1, :, :] = 0
    R[:, 0, :] = 0
    R[:, -1, :] = 0
    R[:, :, 0] = 0
    R[:, :, -1] = 0
    return R.ravel()


sol = newton_krylov(func, s.phi.ravel().copy(), iter=60, verbose=False,
                    f_tol=1e-4, f_rtol=1e-10)
s.phi[:] = sol.reshape(s.nx, s.ny, s.nz)
k0 = s.k0
ite = s.wing_ix[-1]
ku1 = min(k0 + 1, s.nz - 1)
ku2 = min(k0 + 2, s.nz - 1)
kl1 = max(k0 - 1, 0)
kl2 = max(k0 - 2, 0)
zp1, zp2 = s.z[ku1], s.z[ku2]
zm1, zm2 = s.z[kl1], s.z[kl2]
wu = abs(zp2) / max(abs(zp2 - zp1), 1e-12)
wl = abs(zm2) / max(abs(zm2 - zm1), 1e-12)
jmid = len(s.wing_jy) // 2
jm = s.wing_jy[jmid]
print('prescribed G[mid]=%.4f' % s.wake_G[jm])
print('phi above TE col:', np.round(s.phi[ite, jm, ku1], 4),
      np.round(s.phi[ite, jm, ku2], 4), 'z:', round(zp1, 4), round(zp2, 4))
print('phi below TE col:', np.round(s.phi[ite, jm, kl1], 4),
      np.round(s.phi[ite, jm, kl2], 4), 'z:', round(zm1, 4), round(zm2, 4))
pup = wu * s.phi[ite, jm, ku1] - (wu - 1) * s.phi[ite, jm, ku2]
plo = wl * s.phi[ite, jm, kl1] - (wl - 1) * s.phi[ite, jm, kl2]
print('extrapolated jump=%.4f' % (pup - plo))
print('phi at cut faces k0-1,k0:', np.round(s.phi[ite, jm, k0 - 1], 4),
      np.round(s.phi[ite, jm, k0], 4))
