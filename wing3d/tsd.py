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
    def core(n, lo, hi, f0, f1, frac=0.6, le_cluster=True,
             center_cluster=False):
        nc = max(8, int(n * frac))
        u = np.linspace(0, 1, nc)
        if center_cluster:
            # dense at core CENTER (wing/wake plane), sparse at ends:
            # sine remap, slope (1+A*cos(2*pi*u)): A>0 dips at center
            A = 0.75
            u = u + A * np.sin(2 * np.pi * u) / (2 * np.pi)
        elif le_cluster:
            # cosine clustering toward both ends of core (LE/TE emphasis)
            u = 0.5 * (1 - np.cos(np.pi * u))
            u = (u - u.min()) / max(u.max() - u.min(), 1e-12)
        core = f0 + (f1 - f0) * u
        nl = max(2, (n - nc) // 2)
        nr = max(2, n - nc - nl)
        # NOTE (fix): left part must INCREASE lo -> f0 (was assembled
        # backward, folding the grid with negative dx and blowing up SOR)
        left = lo + (f0 - lo) * np.geomspace(0.05, 1, nl) if nl > 1 else []
        right = f1 + (hi - f1) * np.geomspace(0.05, 1, nr) if nr > 1 else []
        # deduplicate junctions (left ends at f0 = core start -> dx=0 ->
        # NaN) and pin exact endpoints before resampling
        g = np.unique(np.concatenate(
            [[lo], np.atleast_1d(left), core, np.atleast_1d(right), [hi]]))
        # resample to exactly n points (monotone interp)
        return np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(g)), g)
    x = core(nx, xlim[0], xlim[1], -0.5 * chord, 2.0 * chord, le_cluster=True)
    y = core(ny, ylim[0], ylim[1], -span / 2 - 0.5, span / 2 + 0.5,
             le_cluster=True)  # cluster at core ends = wing tips
    z = core(nz, zlim[0], zlim[1], -0.5 * chord, 0.5 * chord,
             le_cluster=False, center_cluster=True)
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
                 farfield=None, wake_sign=+1.0, linear=False,
                 bound_ramp=True):
        self.Minf = Minf
        self.linear = linear  # freeze A = 1-Min^2 (diagnostic linear regime)
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
        # Bound-vortex sheet: the wake jump alone realizes only half the
        # lift (no wing-bound circulation). Prescribe the wing jump too,
        # ramped 0 at LE -> Gamma at TE (panel-mu-exact version is the
        # refinement; linear ramp already recovers the missing half).
        self.bound_ramp = bound_ramp and (wake_gamma is not None)
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
        """Conservative TSD residual, compact finite-volume form.

        Node divergence straight from face fluxes (no face-averaging:
        averaging widens the stencil to i+-2 with a checkerboard null
        space that blows up SOR). Murman-Cole switching picks the
        upwind flux at supersonic faces.
        """
        phi, nx, ny, nz = self.phi, self.nx, self.ny, self.nz
        dx, dy, dz = self.dx, self.dy, self.dz
        R = np.zeros_like(phi)
        # ---- x ----
        phix_f = (phi[1:, :, :] - phi[:-1, :, :]) / dx[:, None, None]
        if self.linear:
            A_f = np.full_like(phix_f, max(1 - self.Minf ** 2, 0.05))
            Fx = A_f * phix_f
        else:
            A_f = self._coeff_A(phix_f)
            Fx = A_f * phix_f
            # Murman-Cole: backward flux at supersonic faces
            sub = A_f < 0
            Fup = np.zeros_like(Fx)
            Fup[1:] = A_f[1:] * (phi[1:-1, :, :] - phi[:-2, :, :]) / \
                dx[:-1, None, None]
            Fx = np.where(sub, Fup, Fx)
        # control-volume widths (node-centered)
        wx = np.empty(nx)
        wx[1:-1] = 0.5 * (dx[:-1] + dx[1:])
        wx[0] = dx[0]
        wx[-1] = dx[-1]
        # faces 0..nx-2; node i uses faces i-1, i
        R[1:-1, :, :] += (Fx[1:, :, :] - Fx[:-1, :, :]) / \
            wx[1:-1, None, None]
        # ---- y (central) ----
        phiy_f = (phi[:, 1:, :] - phi[:, :-1, :]) / dy[None, :, None]
        wy = np.empty(ny)
        wy[1:-1] = 0.5 * (dy[:-1] + dy[1:])
        wy[0] = dy[0]
        wy[-1] = dy[-1]
        R[:, 1:-1, :] += (phiy_f[:, 1:, :] - phiy_f[:, :-1, :]) / \
            wy[None, 1:-1, None]
        # ---- z (central) + wake branch cut ----
        phiz_f = (phi[:, :, 1:] - phi[:, :, :-1]) / dz[None, None, :]
        k0 = self.k0
        if len(self.wake_ix) and len(self.wake_jy) and k0 >= 1:
            WI, WJ = np.meshgrid(self.wake_ix, self.wake_jy, indexing='ij')
            phiz_f[WI, WJ, k0 - 1] -= (self.wake_sign * self.wake_scale
                                       * self.wake_G[WJ] / self.dz[k0 - 1])
        if self.bound_ramp and len(self.wing_ix) and len(self.wing_jy) \
                and k0 >= 1:
            BI, BJ = np.meshgrid(self.wing_ix, self.wing_jy, indexing='ij')
            xc = np.clip(self.xc[self.wing_ix] / self.chord, 0.0, 1.0)
            ramp = xc[:, None]  # 0 at LE -> 1 at TE
            phiz_f[BI, BJ, k0 - 1] -= (self.wake_sign * self.wake_scale
                                       * ramp * self.wake_G[BJ]
                                       / self.dz[k0 - 1])
        wz = np.empty(nz)
        wz[1:-1] = 0.5 * (dz[:-1] + dz[1:])
        wz[0] = dz[0]
        wz[-1] = dz[-1]
        R[:, :, 1:-1] += (phiz_f[:, :, 1:] - phiz_f[:, :, :-1]) / \
            wz[None, None, 1:-1]
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
        # exact compact-stencil diagonal |A_ii| (SOR needs the true
        # diagonal for the omega<2 guarantee on stretched grids)
        dx, dy, dz = self.dx, self.dy, self.dz
        wx = np.empty(nx)
        wx[1:-1] = 0.5 * (dx[:-1] + dx[1:])
        wx[0] = dx[0]
        wx[-1] = dx[-1]
        wy = np.empty(ny)
        wy[1:-1] = 0.5 * (dy[:-1] + dy[1:])
        wy[0] = dy[0]
        wy[-1] = dy[-1]
        wz = np.empty(nz)
        wz[1:-1] = 0.5 * (dz[:-1] + dz[1:])
        wz[0] = dz[0]
        wz[-1] = dz[-1]
        beta2 = max(1 - self.Minf ** 2, 0.05)
        ax = np.empty(nx)
        ax[1:-1] = beta2 * (1 / (dx[:-1] * wx[1:-1]) + 1 / (dx[1:] * wx[1:-1]))
        ax[0] = beta2 * 2 / dx[0] ** 2
        ax[-1] = beta2 * 2 / dx[-1] ** 2
        ay = np.empty(ny)
        ay[1:-1] = (1 / (dy[:-1] * wy[1:-1]) + 1 / (dy[1:] * wy[1:-1]))
        ay[0] = 2 / dy[0] ** 2
        ay[-1] = 2 / dy[-1] ** 2
        az = np.empty(nz)
        az[1:-1] = (1 / (dz[:-1] * wz[1:-1]) + 1 / (dz[1:] * wz[1:-1]))
        az[0] = 2 / dz[0] ** 2
        az[-1] = 2 / dz[-1] ** 2
        D = (ax[:, None, None] + ay[None, :, None] + az[None, None, :])
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
