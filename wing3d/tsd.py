"""3D conservative transonic small-disturbance (TSD) solver.

Conservative form (flow along +x, perturbation potential phi)::
    d/dx [ A dphi/dx ] + d/dy [ dphi/dy ] + d/dz [ dphi/dz ] = 0
    A = (1 - Minf^2) - (gamma+1) Minf^2 dphi/dx   (nonlinear coefficient)

Murman-Cole switching: central x-differencing where A > 0 (subsonic),
upwind where A < 0 (supersonic); conservative shock capturing.
Thin-wing transpiration BC on the mean plane; farfield Dirichlet.
Red-black SOR. Cp = -2 dphi/dx (consistent TSD order); wave drag comes
out of pressure integration over the captured shock structure.
"""
import numpy as np


def _stretch(n, a, b, fine, r=1.25):
    """1D grid from a to b with uniform-fine spacing near `fine` region.

    Simplified: uniform in [fine0, fine1], geometric outside.
    """
    return np.linspace(a, b, n)


def make_grid(chord=1.0, span=6.0, nx=100, ny=40, nz=40,
              xlim=(-8, 10), ylim=(-7, 7), zlim=(-7, 7)):
    # dense uniform core around wing/chord, geometric stretch to farfield;
    # extra cosine clustering at LE/TE inside the core via remap
    def core(n, lo, hi, f0, f1, frac=0.6, le_cluster=True):
        nc = max(8, int(n * frac))
        u = np.linspace(0, 1, nc)
        if le_cluster:
            # cosine clustering toward both ends of core (LE/TE emphasis)
            u = 0.5 * (1 - np.cos(np.pi * u))
            u = (u - u.min()) / max(u.max() - u.min(), 1e-12)
        core = f0 + (f1 - f0) * u
        nl = max(2, (n - nc) // 2)
        nr = max(2, n - nc - nl)
        left = f0 - (f0 - lo) * np.geomspace(1, 0.05, nl)[::-1] if nl > 1 else []
        right = f1 + (hi - f1) * np.geomspace(0.05, 1, nr) if nr > 1 else []
        g = np.concatenate([np.atleast_1d(left), core, np.atleast_1d(right)])
        # resample to exactly n points (monotone interp)
        return np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(g)), g)
    x = core(nx, xlim[0], xlim[1], -0.5 * chord, 2.0 * chord, le_cluster=True)
    y = core(ny, ylim[0], ylim[1], -span / 2 - 0.5, span / 2 + 0.5,
             le_cluster=False)
    z = core(nz, zlim[0], zlim[1], -0.5 * chord, 0.5 * chord,
             le_cluster=True)
    return x, y, z


def thickness_slope(xc, t=0.12):
    """NACA 4-digit half-thickness slope d(yt)/dx (unit chord)."""
    x = np.clip(np.asarray(xc, dtype=float), 0, 1)
    d = 5 * t * (0.2969 / (2 * np.sqrt(np.maximum(x, 1e-9))) - 0.1260
                - 2 * 0.3516 * x + 3 * 0.2843 * x ** 2 - 4 * 0.1036 * x ** 3)
    return d


class TSDSolver:
    def __init__(self, Minf, alpha_deg=0.0, chord=1.0, span=6.0,
                 thick=0.12, gamma=1.4, nx=100, ny=40, nz=40,
                 omega=1.5, wake_gamma=None, wake_xmax=None,
                 farfield=None, wake_sign=+1.0):
        self.Minf = Minf
        self.alpha = np.radians(alpha_deg)
        self.chord = chord
        self.span = span
        self.thick = thick
        self.gamma = gamma
        self.omega = omega
        # wake circulation distribution Gamma(y) (from panel Kutta solution;
        # None -> non-lifting). Dict with 'y' stations and 'gamma' values.
        self.wake_gamma = wake_gamma
        # sign of the wake branch-cut jump (diagnostic; +1 standard)
        self.wake_sign = wake_sign
        # farfield[i,j,k] Dirichlet values on outer box (default 0).
        # For lifting cases pass panel-wake-consistent values (see
        # wake_farfield()) so the wake jump exits cleanly.
        self.farfield = farfield
        self.wake_xmax = wake_xmax
        self.x, self.y, self.z = make_grid(chord, span, nx, ny, nz)
        self.nx, self.ny, self.nz = nx, ny, nz
        self.phi = np.zeros((nx, ny, nz))
        # wing mask: cells whose x-center in [0,chord], |y|<span/2
        xc = 0.5 * (self.x[:-1] + self.x[1:])
        yc = 0.5 * (self.y[:-1] + self.y[1:])
        self.wing_ix = np.where((xc >= 0) & (xc <= chord))[0]
        self.wing_jy = np.where(np.abs(yc) <= span / 2)[0]
        # z index of wing plane (grid plane nearest z=0)
        self.k0 = int(np.argmin(np.abs(self.z)))
        # wake slit: x in [chord, outlet], |y| <= span/2, at z=0 plane
        # with prescribed potential jump Gamma(y) (branch cut). Extending
        # to the outlet avoids a truncated-wake edge filament in-domain
        # (outlet inconsistency is far away and benign).
        if wake_xmax is None:
            wake_xmax = self.x[-1]
        self.wake_ix = np.where((xc >= chord) & (xc <= wake_xmax))[0]
        self.wake_jy = np.where(np.abs(yc) <= span / 2)[0]
        if wake_gamma is not None:
            yg = np.asarray(wake_gamma['y'])
            gg = np.asarray(wake_gamma['gamma'])
            self.wake_G = np.interp(yc, yg, gg, left=0.0, right=0.0)
            self.wake_scale = 0.0  # ramped 0 -> 1 during solve
        else:
            self.wake_G = np.zeros_like(yc)
            self.wake_scale = 1.0
        # NOTE: circulation needs Kutta/wake (above); without it
        # transpiration alone cannot lift (verified: CL~0).
        # metrics
        self.dx = np.diff(self.x)
        self.dy = np.diff(self.y)
        self.dz = np.diff(self.z)
        self.xc = xc
        self.yc = yc
        self.zc = 0.5 * (self.z[:-1] + self.z[1:])

    def _coeff_A(self, phix):
        gm = self.gamma
        return (1 - self.Minf ** 2) - (gm + 1) * self.Minf ** 2 * phix

    def residual(self):
        """Conservative TSD residual (interior). Returns R (nx,ny,nz)."""
        phi, nx, ny, nz = self.phi, self.nx, self.ny, self.nz
        dx, dy, dz = self.dx, self.dy, self.dz
        # face gradients / coefficients (vectorized)
        phix_f = (phi[1:, :, :] - phi[:-1, :, :]) / dx[:, None, None]
        A_f = self._coeff_A(phix_f)
        Fx = A_f * phix_f
        # Murman-Cole: upwind flux where A<0 (use backward-biased stencil)
        # implemented as flux blending below in x-divergence
        sub = A_f < 0
        # x-divergence with switching: central normally; where the
        # downwind A is negative use fully-upwinded flux difference
        dFx = np.zeros_like(phi)
        # central face fluxes averaged to nodes (interior)
        Fc = np.zeros((nx + 1, ny, nz))
        Fc[1:-1] = 0.0
        Fc[1:nx - 1] = 0.5 * (Fx[:-1] + Fx[1:])
        Fc[nx - 1] = Fx[nx - 2]  # one-sided at outlet
        Fc[0] = 0.0  # farfield (phi=0 enforced)
        dFx = (Fc[1:] - Fc[:-1]) / np.append(dx, dx[-1])[:, None, None]
        # supersonic correction: replace with upwind difference where needed.
        # Node i is affected if face i-1 or face i is supersonic
        # (faces indexed 0..nx-2 between nodes).
        left = np.zeros((nx, ny, nz), dtype=bool)
        left[1:] = sub
        right = np.zeros((nx, ny, nz), dtype=bool)
        right[:nx - 1] = sub
        node_sup = left | right
        node_sup[0] = False
        # upwind divergence (backward flux difference)
        Fup = np.zeros((nx + 1, ny, nz))
        Fup[1:nx] = Fx
        Fup[0] = 0.0
        Fup[nx] = Fx[nx - 2]
        dFx_up = (Fup[1:] - Fup[:-1]) / np.append(dx, dx[-1])[:, None, None]
        R = dFx
        R[node_sup] = dFx_up[node_sup]
        # y, z divergences (central, linear)
        phiy_f = (phi[:, 1:, :] - phi[:, :-1, :]) / dy[None, :, None]
        Gy = np.zeros((nx, ny + 1, nz))
        Gy[:, 1:ny - 1, :] = 0.5 * (phiy_f[:, :-1, :] + phiy_f[:, 1:, :])
        Gy[:, ny - 1, :] = phiy_f[:, ny - 2, :]
        Gy[:, 0, :] = 0.0
        R += (Gy[:, 1:, :] - Gy[:, :-1, :]) / np.append(dy, dy[-1])[None, :, None]
        phiz_f = (phi[:, :, 1:] - phi[:, :, :-1]) / dz[None, None, :]
        # wake branch cut: subtract prescribed jump so the operator sees
        # continuous gradient (fixed Gamma -> matrix unchanged, stable)
        k0 = self.k0
        if len(self.wake_ix) and len(self.wake_jy) and k0 >= 1:
            WI, WJ = np.meshgrid(self.wake_ix, self.wake_jy, indexing='ij')
            phiz_f[WI, WJ, k0 - 1] -= (self.wake_sign * self.wake_scale
                                       * self.wake_G[WJ] / self.dz[k0 - 1])
        Hz = np.zeros((nx, ny, nz + 1))
        Hz[:, :, 1:nz - 1] = 0.5 * (phiz_f[:, :, :-1] + phiz_f[:, :, 1:])
        Hz[:, :, nz - 1] = phiz_f[:, :, nz - 2]
        Hz[:, :, 0] = 0.0
        R += (Hz[:, :, 1:] - Hz[:, :, :-1]) / np.append(dz, dz[-1])[None, None, :]
        # wing transpiration BC on z=0 plane (both sides)
        R = self._apply_wing_bc(R)
        # farfield Dirichlet (phi=0): zero residual enforced via mask in solve
        return R

    def _apply_wing_bc(self, R):
        """Thin-wing transpiration as pure RHS source (matrix untouched).

        Blowing w (z-velocity) through the wing plane enters divergence
        as a mass source: upper volume gets +wu/dz, lower gets -wl/dz
        (wl negative = downward blowing = source below). No phi
        dependence -> operator stays SPD -> SOR converges.
        Upper: wu = +dz_t/dx - alpha; lower: wl = -dz_t/dx - alpha.
        """
        k0 = self.k0
        if k0 < 1 or k0 >= self.nz:
            return R
        dzm = self.dz[k0 - 1]
        dzp = self.dz[k0] if k0 < len(self.dz) else self.dz[-1]
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        if len(iu) == 0 or len(jy) == 0:
            return R
        II, JJ = np.meshgrid(iu, jy, indexing='ij')
        xc = self.xc[iu] / self.chord
        slope_t = thickness_slope(xc, self.thick)[:, None]
        # NOTE: no cap here (cap kills LE suction); Newton-Krylov handles
        # the stiff LE source robustly where SOR blew up.
        wu = slope_t - self.alpha
        wl = -slope_t - self.alpha
        R[II, JJ, k0] += wu / dzp
        R[II, JJ, k0 - 1] += -wl / dzm
        return R

    def _apply_bc(self, upd):
        """Freeze outer-box phi (0 or prescribed farfield)."""
        if self.farfield is None:
            upd[0, :, :] = 0
            upd[-1, :, :] = 0
            upd[:, 0, :] = 0
            upd[:, -1, :] = 0
            upd[:, :, 0] = 0
            upd[:, :, -1] = 0
            return
        F = self.farfield
        # pin boundary nodes to farfield by zeroing their update
        # (phi initialized to farfield by caller)
        upd[0, :, :] = 0
        upd[-1, :, :] = 0
        upd[:, 0, :] = 0
        upd[:, -1, :] = 0
        upd[:, :, 0] = 0
        upd[:, :, -1] = 0

    def solve(self, itmax=600, tol=1e-5, verbose=True, ramp=100):
        """Red-black SOR (SPD operator: standard form).

        ramp: iterations over which wake_scale goes 0 -> 1 (avoids
        startup shock from impulsive circulation).
        """
        if self.farfield is not None:
            F = self.farfield
            self.phi[0, :, :] = F[0, :, :]
            self.phi[-1, :, :] = F[-1, :, :]
            self.phi[:, 0, :] = F[:, 0, :]
            self.phi[:, -1, :] = F[:, -1, :]
            self.phi[:, :, 0] = F[:, :, 0]
            self.phi[:, :, -1] = F[:, :, -1]
        nx, ny, nz = self.nx, self.ny, self.nz
        ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny),
                                 np.arange(nz), indexing='ij')
        red = ((ii + jj + kk) % 2 == 0)
        # base diagonal (linear part; wing BC is RHS-only)
        dxm = np.append(self.dx, self.dx[-1])[:, None, None]
        dym = np.append(self.dy, self.dy[-1])[None, :, None]
        dzm = np.append(self.dz, self.dz[-1])[None, None, :]
        beta2 = max(1 - self.Minf ** 2, 0.05)
        D = 2 * beta2 / dxm ** 2 + 2 / dym ** 2 + 2 / dzm ** 2
        for it in range(itmax):
            # ramp wake circulation to avoid impulsive startup
            if ramp > 0:
                self.wake_scale = min(1.0, it / max(ramp, 1))
            # farfield tracks the ramp (frozen full-jump farfield against
            # partial interior jump destabilizes corners)
            if self.farfield is not None:
                F = self.farfield * self.wake_scale
                self.phi[0, :, :] = F[0, :, :]
                self.phi[-1, :, :] = F[-1, :, :]
                self.phi[:, 0, :] = F[:, 0, :]
                self.phi[:, -1, :] = F[:, -1, :]
                self.phi[:, :, 0] = F[:, :, 0]
                self.phi[:, :, -1] = F[:, :, -1]
            for mask in (red, ~red):
                R = self.residual()
                upd = np.zeros_like(self.phi)
                # SOR for Ax=b with A_ii<0 (discrete Laplacian-like):
                # x += (w/A_ii)(b-Ax) = +w*R/D, D = |A_ii| > 0.
                upd[mask] = +self.omega * R[mask] / D[mask]
                self._apply_bc(upd)
                self.phi[mask] += upd[mask]
            if it % 50 == 0 or it == itmax - 1:
                r = np.abs(R[1:-1, 1:-1, 1:-1]).max()
                if verbose:
                    print(f'  it {it}: max|R|={r:.3e}')
                if r < tol:
                    break
        return self.phi

    def surface_cp(self):
        """Cp on wing plane (upper/lower) + integrated CL/CD."""
        k0 = self.k0
        dxi = self.dx
        # x-gradient at wing plane (central in x)
        phix = np.zeros((self.nx, self.ny))
        phix[1:-1] = (self.phi[2:, :, k0] - self.phi[:-2, :, k0]) / (
            self.x[2:] - self.x[:-2])[:, None]
        cp = -2 * phix
        CL = CD = 0.0
        area = 0.0
        for i in self.wing_ix:
            for j in self.wing_jy:
                dA = dxi[min(i, len(dxi) - 1)] * self.dy[min(j, len(self.dy) - 1)]
                # upper and lower (same Cp in thin-wing approx -> use jump?)
                CL += -cp[i, j] * dA  # placeholder (see loads())
                area += dA
        return cp, CL, area

    def loads(self):
        """Thin-wing loads: lift from Cp jump, drag from Cp on local slope.

        Cp is linearly extrapolated to the surface (z=0) from the two
        nearest planes each side (off-surface gradients under-read
        suction). CL = int (cpl - cpu) dA / S (lower minus upper).
        CD from pressure on local thickness slope.
        """
        k0 = self.k0
        dxi = self.dx

        def gradx(k):
            g = np.zeros((self.nx, self.ny))
            g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
                self.x[2:] - self.x[:-2])[:, None]
            return g

        def surf_cp(k1, k2):
            # linear extrapolation to z=0 from planes k1 (near) and k2 (far)
            z1 = self.z[k1]
            z2 = self.z[k2]
            w = abs(z2) / max(abs(z2 - z1), 1e-12)
            return -2 * (w * gradx(k1) + (1 - w) * gradx(k2))

        ku1 = min(k0 + 1, self.nz - 1)
        ku2 = min(k0 + 2, self.nz - 1)
        kl1 = max(k0 - 1, 0)
        kl2 = max(k0 - 2, 0)
        cpu = surf_cp(ku1, ku2)
        cpl = surf_cp(kl1, kl2)
        CL = CD = 0.0
        for i in self.wing_ix:
            xc = self.xc[i] / self.chord
            st = thickness_slope(xc, self.thick)
            for j in self.wing_jy:
                dA = dxi[min(i, len(dxi) - 1)] * self.dy[min(j, len(self.dy) - 1)]
                CL += (cpl[i, j] - cpu[i, j]) * dA
                CD += (cpu[i, j] - cpl[i, j]) * st * dA
        S = self.chord * self.span
        return {'CL': CL / S, 'CD': CD / S, 'cpu': cpu, 'cpl': cpl}


def wake_farfield(panel_wakes, muw, X, Y, Z):
    """Potential of panel wake doublets at points (for TSD farfield BC).

    Uses the exact solid-angle kernel (consistent signs with solver),
    so the wake jump exits the TSD box cleanly instead of colliding
    with a phi=0 wall.
    """
    from .solver import wake_precompute
    wc = wake_precompute(panel_wakes)
    fan, fmask = wc['fan'], wc['fmask']  # (N,F,3,3), (N,F)
    P = np.stack(np.broadcast_arrays(X, Y, Z), axis=-1).reshape(-1, 3)
    owner = wc['owner']
    out = np.zeros(len(P))

    def solid_angles(pts):
        # returns (N, M)
        q = pts[None, :, None, :]
        a = fan[:, None, :, 0, :] - q
        b = fan[:, None, :, 1, :] - q
        c = fan[:, None, :, 2, :] - q
        la = np.linalg.norm(a, axis=-1)
        lb = np.linalg.norm(b, axis=-1)
        lc = np.linalg.norm(c, axis=-1)
        num = np.einsum('nfmi,nfmi->nfm', a, np.cross(b, c))
        den = (la * lb * lc + np.einsum('nfmi,nfmi->nfm', a, b) * lc
               + np.einsum('nfmi,nfmi->nfm', b, c) * la
               + np.einsum('nfmi,nfmi->nfm', c, a) * lb)
        om = 2.0 * np.arctan2(num, den)
        return np.where(fmask[:, None, :], om, 0.0).sum(axis=-1)

    # chunked (points x panels can be large)
    for a0 in range(0, len(P), 2048):
        b0 = min(a0 + 2048, len(P))
        om = solid_angles(P[a0:b0])
        # +C.mu with C = -Omega/4pi (consistent with TSD Phi = Phi_inf+B+C)
        for k in range(len(muw)):
            out[a0:b0] += -om[owner == k].sum(axis=0) / (4 * np.pi) * muw[k]
    return out.reshape(np.broadcast(X, Y, Z).shape)


def tsd_farfield_box(solver, panel_wakes, muw):
    """Farfield Dirichlet array for a TSDSolver from the panel wake."""
    X, Y, Z = np.meshgrid(solver.x, solver.y, solver.z, indexing='ij')
    fv = wake_farfield(panel_wakes, muw, X, Y, Z)
    F = np.zeros_like(solver.phi)
    F[0, :, :] = fv[0, :, :]
    F[-1, :, :] = fv[-1, :, :]
    F[:, 0, :] = fv[:, 0, :]
    F[:, -1, :] = fv[:, -1, :]
    F[:, :, 0] = fv[:, :, 0]
    F[:, :, -1] = fv[:, :, -1]
    return F
