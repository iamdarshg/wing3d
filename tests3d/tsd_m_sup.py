"""TSD full-potential at M>=1: sup inflow-pinned, outflow-free BC test.

Hypothesis from RCA: at M=1.0 the TSD stays on the subsonic branch
(CL ~ 0.115 incompressible-like) because ALL box faces are free-floating
early-branch DOFs. At M>=1 the flow is hyperbolic-in-x: the upstream
face must be pinned to uniform (phi=0), the downstream face outflow-free,
so the solution can develop the bow-shock/supersonic-pocket character
(near-symmetric, CL collapses, wave drag grows).

Runs M=1.0 and M=1.1 AR6, reports CL/net section + supersonic fraction.
Truth references (OF AR1.6 anchors, NOT AR-matched): M1.0 CL=0.0787,
M1.1 CL=0.0533. AR6 should sit ABOVE the AR1.6 anchor (tip relief).
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_coupled import build_case
from scipy.optimize import newton_krylov


def run_sup(M, span=6.0, epochs=10, free_k=True, profile=None, section=None,
            base=None):
    if base is None:
        s, res, mesh = build_case(M, 1.25, False, True, None, 80, 32, 32,
                                  span, section=section)
    else:
        # reuse a prebuilt state: panel/imprint is Mach-independent, only
        # the operator coefficients depend on Minf (read live in residual)
        s = base
    s.Minf = M
    if profile is not None:
        s.profile = profile
    # strip panel imprints that lock the subsonic circulation:
    # keep wing_jump (ramp shape) since free-k proved unstable, but
    # force new BCs
    s.farfield = None
    s.phi[:] = 0.0
    s.k[:] = 1.0
    nx, ny, nz = s.nx, s.ny, s.nz
    nphi = (nx - 1) * ny * nz   # interior + outflow face (inflow excluded)
    nk = len(s.k)

    def unpack(y):
        phi = np.zeros((nx, ny, nz))
        phi[1:, :, :] = y[:nphi].reshape(nx - 1, ny, nz)
        # inflow face: uniform supersonic freestream (phi=0)
        phi[0, :, :] = 0.0
        return phi

    def func(y):
        s.phi[:] = unpack(y)
        if free_k:
            s.k[:] = y[nphi:]
        R = s.residual()
        # out-flow face free (zero residual row), lateral/top/bottom free
        R[-1, :, :] = 0
        R[:, 0, :] = 0
        R[:, -1, :] = 0
        R[:, :, 0] = 0
        R[:, :, -1] = 0
        R[0, :, :] = 0   # inflow pinned by unpack, no equation needed
        Rk = s.kutta_residual(upstream=2) * 10.0
        # only return the free DOF rows (inflow plane excluded)
        if free_k:
            return np.concatenate([R[1:].ravel(), Rk])
        return R[1:].ravel()

    for ep in range(epochs):
        s._frozen = s._sup_mask()
        s._frozenA = s._face_coeff()
        dx = s.dx
        phix = (s.phi[1:, :, :] - s.phi[:-1, :, :]) / dx[:, None, None]
        s._frozenRho = s._face_rho(phix)
        y0 = np.concatenate([s.phi[1:].ravel(), s.k.copy()]) if free_k \
            else s.phi[1:].ravel()
        y = newton_krylov(func, y0, iter=80, verbose=False,
                          f_tol=1e-4, f_rtol=1e-7)
        s.phi[:] = unpack(y[:nphi]) if free_k else unpack(y)
        if free_k:
            s.k[:] = y[nphi:]
        L = s.fp_loads()
        supn = int(s._sup_mask().sum())
        print('  ep %d: CL=%.4f k_mid=%.3f nsup=%d' % (
            ep, L['CL'], s.k[len(s.k) // 2], supn), flush=True)
    s._frozen = s._frozenA = s._frozenRho = None
    L = s.fp_loads()
    jy = np.asarray(s.wing_jy)
    j0 = jy[np.argmin(np.abs(s.yc[jy]))]
    ix = np.asarray(s.wing_ix)
    xc = s.xc[ix] / s.chord
    cl = float(np.trapezoid((L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc))
    print('M=%.2f: full-wing CL=%.4f  centerline cl=%.4f  Cfmax=%.2f' % (
        M, L['CL'], cl, L['cpu'][ix, j0].max()), flush=True)
    return L, xc, ix, j0


if __name__ == '__main__':
    for M in [1.0, 1.1]:
        print('=== M=%.2f (AR6, inflow-pinned/outflow-free) ===' % M,
              flush=True)
        run_sup(M)