"""Fast epoch solve: 36-colour FD Jacobian + sparse GMRES.

Replaces the FD Newton-Krylov inner solve of ``tsd_m_sup.run_sup``, which
spends ~2.5-5k residual evaluations per epoch and still stops at its
81-step cap with |F| ~ 1e-3 (measured; perf audit 2026-10-22).

Why this works
--------------
The frozen epoch operator is affine in phi (all coefficients lagged) and
its stencil is contained in the box

    di in {-2, -1, 0, +1},  dj in {-1, 0, +1},  dk in {-1, 0, +1}

(the -2 comes from the upwinded x-face at i-1, which reads phi[i-2]).
A product colouring (i mod 4, j mod 3, k mod 3) therefore has 36 colours
such that every row's stencil is single-coloured: 36 residual calls
assemble the exact Jacobian.  On top of that a sparse GMRES with Jacobi
preconditioning solves the epoch to 1e-8 in ~1 s, instead of ~30-340 s.

Boundary rows that ``func`` leaves empty (it zeroes them) are patched
with real BCs instead of being left rank-deficient:

  * lateral / top / bottom box faces: Dirichlet phi = 0 (farfield);
  * outflow plane i = nx-1: linear extrapolation phi[-1] = phi[-2].

This is a deliberate physics change relative to the campaign path (which
left those DOFs to the Krylov trajectory); it is reported separately.

Unknown layout (same as tsd_m_sup.run_sup):
    y = phi[1:, :, :].ravel()   (inflow plane i=0 is pinned to 0)
"""
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import gmres, LinearOperator


def _grid_arrays(nx, ny, nz):
    ii = np.arange(1, nx)[:, None, None]
    jj = np.arange(ny)[None, :, None]
    kk = np.arange(nz)[None, None, :]
    shape = (nx - 1, ny, nz)
    I = np.broadcast_to(ii, shape).copy()
    J = np.broadcast_to(jj, shape).copy()
    K = np.broadcast_to(kk, shape).copy()
    rows = np.arange(I.size).reshape(shape)
    return I, J, K, rows


def _delta(a, idx, period, offsets):
    """offset d such that (idx + d) % period == a, chosen from offsets."""
    m = (a - idx) % period
    out = np.empty_like(m)
    for r, d in enumerate(offsets):
        out[m == r] = d
    return out


def fd_jacobian(func, y0, shape, eps_rel=1e-7, verbose=False):
    """Exact sparse Jacobian of the frozen operator via 36-colour FD.

    Returns (A_csr, F0) where F0 = func(y0).  Rows that func zeroes are
    structurally empty; they are patched by ``patch_bc``.
    """
    nx, ny, nz = shape
    n = y0.size
    I, J, K, rows = _grid_arrays(nx, ny, nz)
    scale = max(1.0, float(np.max(np.abs(y0))))
    eps = eps_rel * scale
    F0 = func(y0)

    rr, cc, vv = [], [], []
    for a in range(4):
        di = _delta(a, I, 4, [0, 1, -2, -1])
        i2 = I + di
        for b in range(3):
            dj = _delta(b, J, 3, [0, 1, -1])
            j2 = J + dj
            for c in range(3):
                dk = _delta(c, K, 3, [0, 1, -1])
                k2 = K + dk
                valid = ((i2 >= 1) & (i2 <= nx - 1) &
                         (j2 >= 0) & (j2 <= ny - 1) &
                         (k2 >= 0) & (k2 <= nz - 1)).ravel()
                if not np.any(valid):
                    continue
                block = rows[(I % 4 == a) & (J % 3 == b) & (K % 3 == c)]
                if block.size == 0:
                    continue
                yp = y0.copy()
                yp[block] += eps
                Fp = func(yp)
                dF = (Fp - F0) / eps
                # dependencies are exact zeros outside the stencil
                m = valid & (dF != 0.0)
                if not np.any(m):
                    continue
                col = ((i2 - 1) * ny + j2) * nz + k2
                rr.append(rows.ravel()[m])
                cc.append(col.ravel()[m])
                vv.append(dF[m])
                if verbose:
                    print('  colour %d%d%d: %d entries' % (
                        a, b, c, int(m.sum())), flush=True)
    A = coo_matrix((np.concatenate(vv),
                    (np.concatenate(rr), np.concatenate(cc))),
                   shape=(n, n)).tocsr()
    A.sum_duplicates()
    return A, F0


def patch_bc(A, b, y0, shape, outflow='extrapolate'):
    """Add the missing boundary equations (rows func leaves empty).

    Dirichlet phi = 0 on j/k box faces; outflow plane linear
    extrapolation.  b is the RHS for the Newton increment Delta = y - y0:
    for a Dirichlet row (phi_new = 0)  ->  b = -y0;
    for extrapolation (phi_new[r] = phi_new[p]) ->
        b = -(y0[r] - y0[p]).
    """
    nx, ny, nz = shape
    n = y0.size
    I, J, K, rows = _grid_arrays(nx, ny, nz)
    rows = rows.ravel()
    lobnd = ((J == 0) | (J == ny - 1) | (K == 0) |
             (K == nz - 1)).ravel()
    outflow_rows = ((I == nx - 1) & ~((J == 0) | (J == ny - 1) |
                                      (K == 0) | (K == nz - 1))).ravel()
    nz_row = A.getnnz(axis=1)
    missing = nz_row == 0
    # any row that is structurally empty and not covered by the two
    # policies (shouldn't happen) falls back to Dirichlet 0
    dirichlet = lobnd | (missing & ~outflow_rows)
    rr, cc, vv = [], [], []
    rd = rows[dirichlet]
    rr.append(rd)
    cc.append(rd)
    vv.append(np.ones(rd.size))
    b[rd] = -y0[rd]
    if outflow == 'extrapolate':
        ro = rows[outflow_rows]
        if ro.size:
            rr.append(ro)
            cc.append(ro)
            vv.append(np.ones(ro.size))
            prev = ro - ny * nz       # (i-1, j, k) in the flat layout
            rr.append(ro)
            cc.append(prev)
            vv.append(-np.ones(ro.size))
            b[ro] = -(y0[ro] - y0[prev])
    elif outflow == 'dirichlet':
        ro = rows[outflow_rows]
        rr.append(ro)
        cc.append(ro)
        vv.append(np.ones(ro.size))
        b[ro] = -y0[ro]
    pat = coo_matrix((np.concatenate(vv),
                      (np.concatenate(rr), np.concatenate(cc))),
                     shape=(n, n)).tocsr()
    return (A + pat).tocsr()


class LinePrecond(object):
    """Exact tridiagonal x-line solve as a GMRES preconditioner.

    The grid anisotropy is 500-1200:1 in x, so Jacobi leaves ~660 Krylov
    iterations.  The x-lines are the dominant coupling: extract the
    tridiagonal part of A along x (diagonal + first off-diagonals,
    dropping the i-2 upwind entry) and solve all (j,k) lines exactly with
    a batched Thomas sweep.  Costs ~1-3 ms/apply.
    """

    def __init__(self, A, shape):
        nx, ny, nz = shape
        n = A.shape[0]
        stride = ny * nz
        Ac = A.tocoo()
        diag = np.zeros(n)
        sub = np.zeros(n)
        sup = np.zeros(n)
        md = Ac.col == Ac.row
        diag[Ac.row[md]] = Ac.data[md]
        ms = Ac.col == Ac.row - stride
        sub[Ac.row[ms]] = Ac.data[ms]
        mu = Ac.col == Ac.row + stride
        sup[Ac.row[mu]] = Ac.data[mu]
        sh = (nx - 1, ny, nz)

        def lines(a):
            return a.reshape(sh).transpose(1, 2, 0).reshape(-1, nx - 1)

        self.shape = shape
        self.n = n
        self.d = lines(diag)
        self.dl = lines(sub)
        self.du = lines(sup)

    def apply(self, v):
        nx, ny, nz = self.shape
        L = nx - 1
        V = v.reshape(nx - 1, ny, nz).transpose(1, 2, 0).reshape(
            -1, L).astype(float, copy=True)
        d, dl, du = self.d, self.dl, self.du
        cp = np.empty_like(V)
        dp = np.empty_like(V)
        den = np.where(np.abs(d[:, 0]) < 1e-300, 1.0, d[:, 0])
        cp[:, 0] = du[:, 0] / den
        dp[:, 0] = V[:, 0] / den
        for i in range(1, L):
            den = d[:, i] - dl[:, i] * cp[:, i - 1]
            den = np.where(np.abs(den) < 1e-300, 1.0, den)
            cp[:, i] = du[:, i] / den
            dp[:, i] = (V[:, i] - dl[:, i] * dp[:, i - 1]) / den
        X = np.empty_like(V)
        X[:, L - 1] = dp[:, L - 1]
        for i in range(L - 2, -1, -1):
            X[:, i] = dp[:, i] - cp[:, i] * X[:, i + 1]
        return X.reshape(ny, nz, L).transpose(2, 0, 1).ravel()


def solve_epoch(func, y0, shape, rtol=1e-8, restart=60, maxiter=1500,
                newton_max=2, outflow='extrapolate', precond='line',
                verbose=False):
    """Solve one frozen epoch operator.

    Returns (y, info) with info = {'maxF0','maxF','gmres_info','n_evals',
    'linear_iters'}.
    """
    n_evals = 0
    y = np.array(y0, dtype=float, copy=True)
    info = {'maxF0': None, 'maxF': None, 'gmres_info': [],
            'n_evals': 0, 'linear_iters': []}
    for it in range(newton_max):
        A, F0 = fd_jacobian(func, y, shape)
        n_evals += 37
        maxF0 = float(np.max(np.abs(F0)))
        if info['maxF0'] is None:
            info['maxF0'] = maxF0
        if verbose:
            print('  newton %d: |F|=%.3e  nnz=%d' % (it, maxF0, A.nnz),
                  flush=True)
        if maxF0 < 1e-10:
            break
        b = -F0.copy()
        A = patch_bc(A, b, y, shape, outflow=outflow)
        M = None
        if precond == 'line':
            try:
                M = LinearOperator(A.shape,
                                   matvec=LinePrecond(A, shape).apply)
            except Exception:
                M = None
        if M is None:
            diag = A.diagonal()
            diag[np.abs(diag) < 1e-300] = 1.0
            M = LinearOperator(A.shape, matvec=lambda v: v / diag)
        n = y.size
        bnorm = np.linalg.norm(b)
        if bnorm == 0.0:
            break
        count = [0]

        def _cb(_r):
            count[0] += 1
            if verbose and count[0] % 200 == 0:
                print('    gmres iters %d (%.0fs)' % (
                    count[0], time.time() - tg0), flush=True)

        tg = time.time()
        tg0 = tg
        delta, ginfo = gmres(A, b / bnorm, rtol=rtol, atol=0.0,
                             restart=restart, maxiter=maxiter, M=M,
                             callback=_cb, callback_type='pr_norm')
        delta = delta * bnorm
        lin_res = float(np.linalg.norm(A @ delta - b) / bnorm)
        info['gmres_info'].append(int(ginfo))
        info['linear_iters'].append(int(count[0]))
        y = y + delta
        F1 = func(y)
        n_evals += 1
        maxF = float(np.max(np.abs(F1)))
        info['maxF'] = maxF
        if verbose:
            print('  newton %d -> |F|=%.3e  gmres info=%d iters=%d '
                  'linres=%.1e  %.1fs' % (
                      it, maxF, ginfo, count[0], lin_res,
                      time.time() - tg), flush=True)
        if maxF < rtol * max(1.0, bnorm):
            break
    info['n_evals'] = n_evals
    return y, info


def _selftest(nx=24, ny=8, nz=8, M=0.85):
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from tsd_coupled import build_case

    for linear, fullpot in [(True, False), (False, True)]:
        tag = 'linear' if linear else 'fullpot'
        s, res, mesh = build_case(M, 1.25, linear, fullpot, None,
                                  nx, ny, nz, 1.6)
        s.farfield = None
        s.phi[:] = 0.0
        s.k[:] = 1.0
        nxi, nyi, nzi = s.nx, s.ny, s.nz

        def unpack(y):
            phi = np.zeros((nxi, nyi, nzi))
            phi[1:] = y.reshape(nxi - 1, nyi, nzi)
            return phi

        def func(y):
            s.phi[:] = unpack(y)
            R = s.residual()
            R[-1] = 0
            R[:, 0] = 0
            R[:, -1] = 0
            R[:, :, 0] = 0
            R[:, :, -1] = 0
            R[0] = 0
            return R[1:].ravel()

        if fullpot:
            s._frozen = s._sup_mask()
            s._frozenA = s._face_coeff()
            phix = (s.phi[1:] - s.phi[:-1]) / s.dx[:, None, None]
            s._frozenRho = s._face_rho(phix)
            s._frozenRhon = s._nodal_rho()
        y0 = s.phi[1:].ravel().copy()
        shape = (nxi, nyi, nzi)
        n = y0.size
        rng = np.random.default_rng(0)

        yt = 0.01 * rng.standard_normal(n)
        F0 = func(np.zeros(n))
        F1 = func(yt)
        F2 = func(2.0 * yt)
        aff = float(np.max(np.abs(F2 - (2.0 * F1 - F0))) /
                    max(np.max(np.abs(F1)), 1e-30))

        A, Fj = fd_jacobian(func, y0, shape)
        v = rng.standard_normal(n)
        t = 1e-6
        d = (func(y0 + t * v) - Fj) / t
        jd = float(np.linalg.norm(d - A @ v) /
                   max(np.linalg.norm(A @ v), 1e-30))

        y, info = solve_epoch(func, y0, shape, verbose=True)
        print('[%s] n=%d nnz=%d  affinity=%.2e  jac_err=%.2e  '
              '|F| %.2e -> %.2e  gmres=%s' % (
                  tag, n, A.nnz, aff, jd, info['maxF0'],
                  info['maxF'], info['gmres_info']), flush=True)


if __name__ == '__main__':
    _selftest()
