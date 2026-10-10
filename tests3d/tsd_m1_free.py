"""Free-circulation TSD at M1.0: no panel imprints.

Drops wing_jump imprint (-> ramp x wake shape, LEVEL set by free k),
panel farfield (-> uniform 0), tare (-> true Kutta: zero TE jump).
Newton-Krylov over [phi, k] must FIND the sonic circulation.
Compares vs OF anchor 0.0787 and imprinted 0.1145.
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_coupled import build_case
from scipy.optimize import newton_krylov

s, res, mesh = build_case(1.0, 1.25, False, True)
# ---- strip the imprints ----
s.wing_jump = None          # ramp x wake_G shape; level via free k
s.farfield = None           # uniform freestream Dirichlet
s.phi[:] = 0.0              # cold start (was panel-box init)
s.k[:] = 1.0
nphi = s.nx * s.ny * s.nz


def func(x):
    s.phi[:] = x[:nphi].reshape(s.nx, s.ny, s.nz)
    s.k[:] = x[nphi:]
    R = s.residual()
    R[0, :, :] = 0
    R[-1, :, :] = 0
    R[:, 0, :] = 0
    R[:, -1, :] = 0
    R[:, :, 0] = 0
    R[:, :, -1] = 0
    Rk = s.kutta_residual(upstream=2) * 10.0  # TRUE Kutta: TE jump = 0
    return np.concatenate([R.ravel(), Rk])


def section_cl(s):
    L = s.fp_loads()
    jy = np.asarray(s.wing_jy)
    j0 = jy[np.argmin(np.abs(s.yc[jy]))]
    ix = np.asarray(s.wing_ix)
    xc = s.xc[ix] / s.chord
    cl = float(np.trapezoid((L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc))
    return L['CL'], cl


for epoch in range(12):
    s._frozen = s._sup_mask()
    s._frozenA = s._face_coeff()
    dx = s.dx
    phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
    s._frozenRho = s._face_rho(phix)
    x0 = np.concatenate([s.phi.ravel(), s.k.copy()])
    sol = newton_krylov(func, x0, iter=80, verbose=False,
                        f_tol=1e-5, f_rtol=1e-8)
    s.phi[:] = sol[:nphi].reshape(s.nx, s.ny, s.nz)
    s.k[:] = sol[nphi:]
    CL, cl = section_cl(s)
    print('  epoch %d: CL=%.4f cl_sec=%.4f k_mid=%.3f nsup=%d' % (
        epoch, CL, cl, s.k[len(s.k) // 2],
        int(s._sup_mask().sum())), flush=True)
s._frozen = s._frozenA = s._frozenRho = None
CL, cl = section_cl(s)
print('FREE: full-wing CL=%.4f centerline cl=%.4f (OF 0.0787)' % (CL, cl))
# ---- field-level diagnosis: pocket structure at centerline ----
gm = s.gamma
jy = np.asarray(s.wing_jy)
j0 = jy[np.argmin(np.abs(s.yc[jy]))]
ix = np.asarray(s.wing_ix)
xc = s.xc[ix] / s.chord
k0 = s.k0
gx_up = (s.phi[ix + 1, j0, k0 + 1] - s.phi[ix - 1, j0, k0 + 1]) / (
    s.x[ix + 1] - s.x[ix - 1])
gx_lo = (s.phi[ix + 1, j0, k0 - 1] - s.phi[ix - 1, j0, k0 - 1]) / (
    s.x[ix + 1] - s.x[ix - 1])


def mach_of_gx(gx):
    q2 = np.maximum((1.0 + gx) ** 2, 1e-9)
    Trat = np.maximum(1 + 0.5 * (gm - 1) * s.Minf ** 2 * (1 - q2), 1e-9)
    return np.sqrt(q2 * s.Minf ** 2 / Trat)


mu, ml = mach_of_gx(gx_up), mach_of_gx(gx_lo)
print(' x   : M-up  M-lo')
for i in range(0, len(ix), 3):
    print(' %.2f : %.2f  %.2f' % (xc[i], mu[i], ml[i]))
np.save('tmp_m1_field.npy',
        np.stack([xc, mu, ml], axis=1))
print('saved field slice')
