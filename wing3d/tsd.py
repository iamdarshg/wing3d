"""3D conservative transonic small-disturbance (TSD) solver.

STATUS: linear regime validated (panel-exact 1.4% via prescribed
panel-mu jumps). Nonlinear (shock) regime PARKED pending bound-vortex
freedom -- design note below.

VORTEX-FREEDOM DESIGN (measured 2026-10-08):
- Panel decomposition (NACA0012 AR6 a=4): body-only CL=0.267,
  wake-only CL=-0.196, full 0.369. The wake UNLOADS the TE (Kutta
  enforcer); it does not lift. Lift lives in the BODY (bound vortex).
- TSD wake-only prescription induces ~half lift (0.203) -- correct
  for a trailed system without its bound mate (Helmholtz), not a bug.
- Prescribing the wing jump too pins the linear state (Kutta then
  auto-satisfied for any scale; tangency monotonic in scale) so
  shocks/lift-growth cannot develop. Kelvin: shocks don't create
  circulation; the Kutta condition reselects Gamma, which must be FREE.
- REQUIRED: explicit bound-vortex sheet unknown gamma_b(x,y) on the
  wing, solved WITH phi: gamma_b from flow-tangency on the wing,
  phi from TSD with gamma_b as jump source, closed by Kutta
  (gamma_b(TE) = wake Gamma, updated). I.e. VLM+TSD coupling.
  Estimated: hours, high value (unlocks transonic wing3d).
- Until then: use OpenFOAM central anchors (cases/openfoam/transonic/)
  for transonic truth; TSD linear + panel jumps for subsonic.

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
import os


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


def flat_plate_jump(xc):
    """Cumulative flat-plate loading, int_0^x A(t) dt.

    A(x) = (2/pi) sqrt((1-x)/x) is the thin-airfoil (flat plate)
    vortex strength, normalised so int_0^1 A dx = 1: the sheet
    amplitude Gamma then IS the total section circulation, and
    dF/dx = A is the classic LE-singular loading. F is the POTENTIAL
    JUMP across the bound sheet (the quantity the branch cut carries),
    so it must be the integral, with F(0)=0 (no jump at the LE) and
    F(1)=1 (the TE jump equals the trailed wake jump: Kutta/Helmholtz).
    Closed form: (2/pi)[sqrt(x(1-x)) + arcsin(sqrt(x))].
    """
    x = np.clip(np.asarray(xc, dtype=float), 0.0, 1.0)
    return (2.0 / np.pi) * (np.sqrt(x * (1.0 - x)) + np.arcsin(np.sqrt(x)))


def section_profile_xy(x, y_up, y_lo, n=400):
    """Build slope/thickness interpolants from an arbitrary airfoil.

    x, y_up, y_lo: chord-normalized coordinates (upper/lower as functions
    of x, LE->TE). Returns dict with 'slope' (d(y_up-y_lo)/2/dx) and
    'half' (half-thickness) callables of x/c, plus the raw points.
    """
    x = np.clip(np.asarray(x, dtype=float), 0, 1)
    y_up = np.asarray(y_up, dtype=float)
    y_lo = np.asarray(y_lo, dtype=float)
    order = np.argsort(x)
    x = x[order]; y_up = y_up[order]; y_lo = y_lo[order]
    # half-thickness (symmetric part) and mean line (camber)
    half = 0.5 * (y_up - y_lo)
    mean = 0.5 * (y_up + y_lo)
    # resample on uniform x/c
    xg = np.linspace(0, 1, n)
    hg = np.interp(xg, x, half)
    mg = np.interp(xg, x, mean)
    slope = np.gradient(hg, xg)
    camber = np.gradient(mg, xg)

    def slope_fn(xc):
        xc = np.clip(np.asarray(xc, dtype=float), 0, 1)
        return np.interp(xc, xg, slope)

    def half_fn(xc):
        xc = np.clip(np.asarray(xc, dtype=float), 0, 1)
        return np.interp(xc, xg, hg)

    def camber_fn(xc):
        """Mean-line (camber line) slope d(yc)/dx."""
        xc = np.clip(np.asarray(xc, dtype=float), 0, 1)
        return np.interp(xc, xg, camber)

    def mean_fn(xc):
        """Mean-line (camber line) height yc(x)."""
        xc = np.clip(np.asarray(xc, dtype=float), 0, 1)
        return np.interp(xc, xg, mg)

    return {'slope': slope_fn, 'half': half_fn, 'camber': camber_fn,
            'mean': mean_fn, 'x': xg, 'half_thickness': hg,
            'mean_line': mg, 'slope_arr': slope, 'camber_arr': camber}


class TSDSolver:
    def __init__(self, Minf, alpha_deg=0.0, chord=1.0, span=6.0,
                 thick=0.12, gamma=1.4, nx=100, ny=40, nz=40,
                 omega=1.5, wake_gamma=None, wake_xmax=None,
                 farfield=None, wake_sign=+1.0, linear=False,
                 bound_ramp=True, wing_jump=None, fullpot=False,
                 box=None, profile=None):
        self.Minf = Minf
        self.linear = linear  # freeze A = 1-Min^2 (diagnostic linear regime)
        # full potential (density-coupled exact potential, valid to M1+)
        # instead of TSD linearization; shares grid/BCs/wake/line-solver
        self.fullpot = fullpot and not linear
        self._frozen = None  # frozen supersonic map (nonlinear epochs)
        self._frozenA = None  # frozen full A coefficient (nonlinear epochs)
        # exact bound jump from panel mu (array over (wing_ix, wing_jy));
        # None -> linear LE->TE ramp of wake Gamma (crude fallback)
        self.wing_jump = wing_jump
        self.alpha = np.radians(alpha_deg)
        self.chord = chord
        self.span = span
        self.thick = thick
        self.profile = profile  # optional section_profile_xy dict
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
        self.x, self.y, self.z = make_grid(chord, span, nx, ny, nz,
                                            **(box or {}))
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
        # FREE bound vortex (Kutta closure): per-span-cell amplitude k
        # scaling the (fixed-shape) wing jump + wake Gamma. Solved with
        # phi via Kutta rows (TE dCp = 0). k=1 reproduces prescribed.
        self.k = np.ones(len(self.wing_jy))
        # ---- FREE bound circulation closed by Kutta (WING3D_KUTTA=1) --
        # gamma_b[j] = per-unit-span circulation at wing span station j,
        # an UNKNOWN of the coupled system. It enters as the potential
        # jump F(x)*gamma_b[j] on the wing (F = cumulative flat-plate
        # loading, 0 at LE -> 1 at TE) and continues unchanged down the
        # wake (Helmholtz), and the lifting camber/alpha transpiration is
        # dropped from _apply_wing_bc so the circulation carries ALL the
        # lift. Closed by the Kutta rows kutta_residual(upstream=2) = 0
        # (TE dCp per span station), added to the residual in
        # tests3d/tsd_m_sup.run_sup.
        # NOTE the name: gamma_b, NOT Gamma -- self.gamma is the ratio
        # of specific heats and a case-collision here would be silent.
        # WING3D_KUTTA=1: transpiration = thickness only, bound sheet
        #   F(x)*Gamma covers the chord, wake trails Gamma  (task spec)
        # WING3D_KUTTA=2: transpiration keeps alpha+camber (the tangency
        #   condition IS the lifting problem), NO bound sheet, the wake
        #   sheet carries the free Gamma (classical TSD/TSFOIL wake).
        # WING3D_KUTTA=3: bound sheet F(x)*Gamma + alpha/camber tangency
        #   both present; the Kutta row selects the sheet amplitude
        #   (the load split is then the sum, as in a VLM+thickness
        #   decomposition).
        # WING3D_KUTTA=4: A10 closure -- bound sheet F(x)*Gamma, thickness-
        #   only transpiration (as mode 1) but the closure row is the
        #   TANGENCY row  wav(x_cp,y_j) - (dyc/dx - alpha) = 0  (the
        #   field's own mean-plane normal velocity, tsd.w_plane, must
        #   balance the geometric slope), plus the MORINO wake condition
        #   Gamma_wake = phi_up(TE) - phi_lo(TE) (measured from the field
        #   instead of prescribed). See tests3d/a10_kutta_report.md.
        # WING3D_KUTTA=5: FREE CHORDWISE DOUBLET DISTRIBUTION (doublet
        #   lattice). A10 measured that one amplitude per span cell cannot
        #   represent a chordwise tangency problem: the selected
        #   circulation was set by the chordwise REDUCTION, not by
        #   physics (pt -0.007 vs gal +0.019 against a cambered LL 0.172).
        #   Here the unknown is the whole potential jump mu(i,j) at nc
        #   chordwise control stations x ng span stations (nc*ng
        #   unknowns), applied as a doublet sheet on the cut. The wake
        #   trails the TE value of mu (Helmholtz). Closure: ONE TANGENCY
        #   ROW PER (chordwise station, span station) --
        #       wav(x_i, y_j) - (dyc/dx(x_i) - alpha) = 0
        #   so every panel feels its OWN camber slope (LE and aft
        #   loading both represented) and the TE row is the Kutta
        #   condition. The chordwise loading shape is then the solution,
        #   not a prescribed F(x). See tests3d/a11_kutta_report.md.
        # WING3D_KUTTA=6: same free chordwise doublet distribution (mu at
        #   nc chordwise control stations x ng span stations) and the same
        #   Helmholtz wake, but the tangency rows are built from the
        #   ANALYTIC DOUBLET-LATTICE KERNEL (tsd.dlm_tangency_residual)
        #   instead of the field-differentiated w_plane estimator that A11
        #   measured to be 3.7x chordwise-biased. D.mu = dyc/dx - alpha is
        #   LINEAR in mu. See tests3d/a12_kutta_report.md.
        # WING3D_KUTTA=7: HYBRID row -- the field-induced normal velocity
        #   with the estimator's chordwise bias subtracted analytically
        #     row = w_field(phi) - [E - D].mu - (dyc/dx - alpha)
        #   E is w_plane's own estimator evaluated on the SAME mu (the
        #   horseshoe field sampled at w_plane's two stencil heights), D is
        #   mode 6's validated kernel. E - D is a FIXED grid/kernel matrix,
        #   so the row is linear in mu and DEPENDS ON phi (mode 6's mu
        #   block was phi-independent -> block-triangular Jacobian -> no
        #   transonic circulation response). See tests3d/a13_kutta_report.md.
        try:
            self.kutta_mode = int(os.environ.get('WING3D_KUTTA', '0') or 0)
        except ValueError:
            self.kutta_mode = 0
        self.kutta_on = self.kutta_mode > 0
        # mode 5/6/7: free chordwise doublet distribution (mu matrix unknown)
        self.kutta_free_mu = self.kutta_mode in (5, 6, 7)
        # mode 6/7: analytic doublet-lattice kernel for the tangency rows
        self.kutta_kernel = self.kutta_mode in (6, 7)
        # mode 7: hybrid field/kernel tangency row (phi-coupled mu block)
        self.kutta_hybrid = self.kutta_mode == 7
        self.kutta_bound = self.kutta_mode in (1, 3, 4, 5, 6, 7)  # sheet from free amps
        self.kutta_thick_only = self.kutta_mode in (1, 4, 5, 6, 7)  # transpiration: dyt
        # closure row family: 1-3 = TE loading dCp (A9), 4 = tangency
        # reduced to one row per span cell, 5/6/7 = tangency per (panel,span)
        self.kutta_tangency = self.kutta_mode in (4, 5, 6, 7)
        # Morino wake condition (mode 4): the wake doublet is the TE
        # potential jump the FIELD carries, not a prescribed copy of
        # gamma_b. WING3D_MORINO=0 turns it back into the prescribed wake
        # (used to measure whether the condition carries information).
        self.kutta_morino = (self.kutta_mode == 4) and (
            os.environ.get('WING3D_MORINO', '1') == '1')
        # mode 5 Kutta row: the last chordwise control station is the TE
        # itself, where the physics is the KUTTA condition (zero loading,
        # dmu/dx -> 0), not tangency on a mean line that ends there.
        # Collocating tangency at the TE instead makes the aft loading
        # satisfy a local condition the continuum does not impose and
        # leaves a TE loading spike (measured: gamma_TE = gamma_max for
        # the dfvlr section). Replacing that row by mu(TE) = mu(last-1)
        # removes the sqrt(1-x) eigensolution -- the standard discrete
        # Kutta in doublet form. WING3D_KUTTA_TE=0 keeps the pure
        # collocation (reported as the ablation).
        self.kutta_te_kutta = (self.kutta_mode in (5, 6, 7)) and (
            os.environ.get('WING3D_KUTTA_TE', '1') == '1')
        # weight of the Kutta row in mu units, commensurate with the
        # tangency rows: MEASURED local self-response d wav / d mu =
        # -6.1 per unit mu (tsd_a11_probe2.py, panels 6 and 11) for the
        # mode-5 field estimator. Mode 6's kernel rows are O(1) in mu, so
        # there the Kutta row is commensurate at weight 1 (it is a pure
        # row scaling either way -- the mu block is linear, so it cannot
        # move the root).
        self.kutta_wte = float(os.environ.get(
            'WING3D_KUTTA_WTE', '1.0' if self.kutta_kernel else '6.1'))
        self.tang_xc = float(os.environ.get('WING3D_TANG_XC', '0.75'))
        # 'gal' = loading-weighted chordwise reduction (default),
        # 'pt'  = single control station at tang_xc
        self.tang_mode = os.environ.get('WING3D_TANG_MODE', 'gal')
        self.gamma_b = np.zeros(len(self.wing_jy))
        # ---- mode 5: free chordwise doublet distribution --------------
        # mu_b[i,j] = potential jump (doublet strength) on chordwise
        # control station i at span station j. Control stations are a
        # SUBSET OF THE FIELD'S OWN chordwise nodes (WING3D_KUTTA_NC,
        # default = all of them, 13 on g80), so the tangency rows read
        # w_plane at nodes instead of interpolating a jump-contaminated
        # estimator. The jump is then mapped back onto every wing node
        # by linear interpolation along the chord.
        self.n_ctrl = 0
        self.i_ctrl = np.zeros(0, dtype=int)
        self.x_ctrl = np.zeros(0)
        self.mu_b = np.zeros((0, len(self.wing_jy)))
        # mode 6 analytic doublet-lattice matrix, built lazily (geometry
        # only -- it does not depend on phi or on mu)
        self._dlmK = None
        self._dlmInfo = None
        # A14: circulation-off REFERENCE field, wav on the wing nodes
        # (shape (n_node, ng)). When set, the hybrid row subtracts its
        # L-mapped downwash: B' = B + w_ref. See set_thickness_reference.
        self._thickRefNodes = None
        # ---- NON-PLANAR BODY SURFACE z_s(x, y) -------------------------
        # The default model is a near-plane: ONE Cartesian z-plane (node
        # k0, sheet on face k0-1) carries the jump + transpiration for the
        # whole span, and the chordwise mask is a single index set shared
        # by every span station. `zsurf` / `wing_mask` / `planform` replace
        # that with a per-(i,j) surface height, a per-(i,j) membership mask
        # and per-span section parameters (chord, LE, twist, section).
        # All of them default to None -> the stock planar path runs and is
        # bit-identical (see `_wing_geom` and tests3d/a16_nonplanar.py).
        self.zsurf = None      # (nx, ny) float: surface height z_s(x_i,y_j)
        self.wing_mask = None  # (nx, ny) bool: cells that carry a body
        self.planform = None   # dict from geometry.*_planform()
        self._geom = None      # cache: snapped per-(i,j) plane geometry
        self.nonplanar = False
        if self.kutta_free_mu:
            self.init_free_mu()
        # (mode 5 mu is initialised from the flat-plate/Lifting-line
        # imprint in tests3d/tsd_m_sup.run_sup, like gamma_b in mode 4)
        # Kutta-row weight, measured against the phi-row Jacobian scale
        # (tests3d/tsd_kutta_weight.py, never an arbitrary 10): the weight
        # that puts |dR_kutt/dGamma| on the phi-row diagonal scale
        # A_scale = 6.33e3 (g80x32x32, M0.3). Measured converged coupled
        # gains |dCp(TE)/dGamma|: 1.18e2 (mode 1, thick-only + bound sheet)
        # -> w=54; 4.72e0 (mode 2, free wake only) -> w=1.5e3.
        # Mode 4 (tangency row): measured converged gain d wav/d Gamma =
        # -0.81 at the 0.75c control point (three converged phi solves,
        # tsd_a10_probe3.py), so w = A_scale/0.81 = 7.24e3/0.81 = 8.9e3.
        # Mode 5 (tangency row from the field-differentiated w_plane) and
        # mode 6 (analytic doublet-lattice kernel): the rows are in
        # velocity units, the phi rows in TSD residual units, so the same
        # A_scale normalisation applies. NOTE the mode-6 mu block is
        # LINEAR and independent of phi, so this weight is pure row
        # conditioning and cannot move the root: kutta_w multiplies both
        # the mu rows and the Kutta row of the same linear block, so it
        # scales rows of a system whose solution is unchanged.
        # Mode 7 (hybrid row) is NOT of that kind: the mu block now depends
        # on phi, so this weight trades the tangency rows against the TSD
        # rows and can move the root. 8.9e3 is inherited from mode 6 (same
        # velocity units, same A_scale construction); the sensitivity of
        # the answer to it is reported as an ablation, not tuned.
        _wd = {1: '50', 2: '1500', 3: '50', 4: '8.9e3', 5: '8.9e3',
               6: '8.9e3', 7: '8.9e3'}.get(self.kutta_mode, '50')
        self.kutta_w = float(os.environ.get('WING3D_KUTTA_W', _wd))

        self._circff = None  # farfield circulation (wake interp target)
        self._frozenRho = None  # frozen face density (fullpot epochs)
        # frozen NODAL density (fullpot epochs). Without this the epoch
        # operator is not affine: rho_n is recomputed from the live phi on
        # every residual call, so "frozen" epochs chase a moving operator
        # (measured affinity drift 1.5e-2; perf audit 2026-10-22).
        self._frozenRhon = None
        # ---- viscous-inviscid interaction (WING3D_BL=1) --------------
        # Integral boundary layer marched on each span strip from the
        # inviscid edge velocity; the displacement thickness delta*
        # feeds back as transpiration (d(delta*)/dx) in _apply_wing_bc.
        # Off by default so the inviscid baseline is bit-identical.
        self.bl_on = os.environ.get('WING3D_BL', '0') == '1'
        self.bl_relax = float(os.environ.get('WING3D_BL_RELAX', '0.5'))
        self.bl_Re = float(os.environ.get('WING3D_RE', '1.8e7'))
        self.bl_xtr = float(os.environ.get('WING3D_BL_XTR', '0.05'))
        self.bl_dstar_u = None  # delta*/c per (wing_ix, wing_jy)
        self.bl_dstar_l = None
        self.bl_hist = []       # one dict per coupling iteration
        self.bl_sep = []        # recorded separation events
        # metrics
        self.dx = np.diff(self.x)
        self.dy = np.diff(self.y)
        self.dz = np.diff(self.z)
        self.xc = xc
        self.yc = yc
        self.zc = 0.5 * (self.z[:-1] + self.z[1:])

    def init_free_mu(self):
        """(Re)build the mode-5 chordwise control-station set.

        Control stations are a subset of the field's own wing nodes
        (WING3D_KUTTA_NC, default = all 13 on g80), so the tangency rows
        read w_plane at nodes rather than interpolating a
        jump-contaminated estimator. mu_b is (n_ctrl, n_span).
        """
        self.kutta_free_mu = True
        self._dlmK = None
        self._dlmInfo = None
        self._dlmE = None      # mode 7: estimator-applied-to-sheet operator
        self._dlmL = None      # mode 7: field-node -> lattice-control map
        self._hybB = None      # mode 7: bias operator L.E - D
        ni = len(self.wing_ix)
        self.n_ctrl = 0
        self.i_ctrl = np.zeros(0, dtype=int)
        self.x_ctrl = np.zeros(0)
        self.mu_b = np.zeros((0, len(self.wing_jy)))
        if ni == 0:
            return
        try:
            nc = int(os.environ.get('WING3D_KUTTA_NC', str(ni)))
        except ValueError:
            nc = ni
        nc = max(2, min(nc, ni))
        idx = np.unique(np.round(np.linspace(0, ni - 1, nc)).astype(int))
        self.i_ctrl = idx
        self.n_ctrl = len(idx)
        # self.xc is built later in __init__, so take the node centres here
        _xc = 0.5 * (self.x[:-1] + self.x[1:])
        self.x_ctrl = _xc[np.asarray(self.wing_ix)[idx]] / self.chord
        self.mu_b = np.zeros((self.n_ctrl, len(self.wing_jy)))

    # ==================================================================
    # NON-PLANAR SURFACE z_s(x, y)   (WING3D_NONPLANAR / set_surface)
    # ==================================================================
    def set_surface(self, zsurf=0.0, mask=None, planform=None, rebuild=True):
        """Install a general body surface z = z_s(x, y) on the wing mask.

        zsurf : (nx, ny) array (node-centred), or a scalar (uniform).
        mask  : (nx, ny) bool of cells carrying a body; None = the stock
                rectangle wing_ix x wing_jy.
        planform : optional dict with per-span section data, see
                `set_planform` (chord / LE / twist / section identity).

        What the generalisation does NOT change: the unknowns (`mu_b`,
        `gamma_b`), the closures, and the grid. What it changes:

          * the sheet (potential jump) and the transpiration are applied at
            the z-FACE that brackets z_s(x_i, y_j) for THAT (i,j), instead
            of at one global face k0-1;
          * the surface quantities (Cp, Cp jump, w_plane, the Kutta rows)
            are extrapolated to z_s(x_i,y_j) instead of to z = 0;
          * the transpiration carries the tangential term of the local
            normal projection, h_y * v (see `_apply_wing_bc`).

        Cartesian snap: z_s is represented by the two adjacent z-node
        stencils with weights (1-theta, theta), theta = the position of z_s
        inside its cell. theta is continuous in z_s, so CL is continuous in
        the dihedral angle; theta = 1 exactly at a node, which is what makes
        the stock z_s = 0 case bit-identical. The residual resolution limit
        is dz near the wing (~2-4% chord on these grids): a height change
        smaller than dz/2 across the semi-span cannot be represented at all.
        """
        import numpy as _np
        if zsurf is None:
            return self.clear_surface()
        self.nonplanar = True
        z = _np.asarray(zsurf, dtype=float)
        if z.ndim == 0:
            z = _np.full((self.nx, self.ny), float(z))
        if z.shape != (self.nx, self.ny):
            raise ValueError('zsurf shape %s != %s' % (z.shape,
                                                       (self.nx, self.ny)))
        self.zsurf = z
        if mask is None:
            self.wing_mask = None
        else:
            m = _np.asarray(mask, dtype=bool)
            if m.shape == (len(self.wing_ix), len(self.wing_jy)):
                full = _np.zeros((self.nx, self.ny), dtype=bool)
                full[_np.ix_(self.wing_ix, self.wing_jy)] = m
                m = full
            if m.shape != (self.nx, self.ny):
                raise ValueError('mask shape %s != %s' % (m.shape,
                                                          (self.nx, self.ny)))
            self.wing_mask = m
        if planform is not None:
            self.planform = planform
        if rebuild:
            self._geom = None
            self._dlmK = self._dlmE = self._dlmL = self._hybB = None
            self._dlmInfo = None
        return self

    def clear_surface(self):
        self.zsurf = None
        self.wing_mask = None
        self.planform = None
        self.nonplanar = False
        self._geom = None
        self._dlmK = self._dlmE = self._dlmL = self._hybB = None
        self._dlmInfo = None
        return self

    def set_planform(self, plan):
        """Rasterise a planform dict (geometry.nonplanar_planform) onto the
        Cartesian (x,y) grid and install it.

        `plan` keys, indexed by span station y:
            y, xle, chord, z, twist_deg, and optionally callables
            slope_t(x_rel, y), slope_c(x_rel, y) for section identity.

        The body membership is per (i,j) CELL: cell (i,j) carries a body
        iff its x-centre lies in [xle(y_j), xle(y_j)+chord(y_j)] for its
        own span station. That single rule is what makes taper, sweep, a
        crank and a centrebody representable -- the stock model applies one
        chordwise index set to every span station, which cannot express any
        of them.
        """
        y = np.asarray(plan['y'], dtype=float)
        xle = np.asarray(plan['xle'], dtype=float)
        ch = np.asarray(plan['chord'], dtype=float)
        zz = np.asarray(plan['z'], dtype=float)
        tw = np.radians(np.asarray(plan.get('twist_deg',
                                            np.zeros_like(y)), dtype=float))
        xc, yc = self.xc, self.yc          # CELL centres, stock convention
        inn_span = np.abs(yc) <= 0.5 * self.span
        xle_j = np.interp(yc, y, xle)
        ch_j = np.interp(yc, y, ch)
        z_cell = np.interp(yc, y, zz)
        tw_j = np.interp(yc, y, tw)
        jset = np.where(inn_span & (ch_j > 0))[0]
        if jset.size == 0:
            raise ValueError('planform has no span stations inside the box')
        lo = float(np.min(xle_j[jset]))
        hi = float(np.max(xle_j[jset] + ch_j[jset]))
        iset = np.where((xc >= lo) & (xc <= hi))[0]
        if iset.size == 0:
            raise ValueError('planform lies outside the grid in x')
        self.wing_ix = iset
        self.wing_jy = jset
        inn = ((xc[iset][:, None] >= xle_j[jset][None, :])
               & (xc[iset][:, None] <= (xle_j[jset] + ch_j[jset])[None, :]))
        z_node = np.interp(self.y, y, zz)
        zsurf = np.repeat(z_node[None, :], self.nx, axis=0)
        meta = {'y': y, 'xle_j': xle_j, 'c_j': ch_j, 'z_j': z_cell,
                'tw_j': tw_j, 'chord': ch, 'xle': xle, 'z': zz}
        for key in ('slope_t', 'slope_c'):
            if key in plan:
                meta[key] = plan[key]
        self.set_surface(zsurf, mask=inn, planform=meta)
        if self.kutta_free_mu:
            self.init_free_mu()
        return self

    def _expand(self, inn):
        """(ni,ng) mask -> full (nx,ny) bool."""
        import numpy as _np
        m = _np.zeros((self.nx, self.ny), dtype=bool)
        m[_np.ix_(self.wing_ix, self.wing_jy)] = inn
        return m

    def _wing_geom(self, rebuild=False):
        """Per-(i,j) snapped-plane geometry (cached).

        Returns dict of (ni, ng) arrays, ii/jj the node index grids:
          ks    node index whose stencil represents the surface
          f     = ks-1: the z-FACE the sheet sits on
          th    in-cell fraction; the stencil of face f is weighted th and
                that of face f-1 is weighted (1-th). th == 1 exactly when
                z_s is a z-node -> the stock scheme, bit for bit.
          zs    surface height, zs_x/zs_y its slopes
          mask  active cells
          w_st  (2, ni, ng) face indices / weights of the two stencils
        """
        import numpy as _np
        if self._geom is not None and not rebuild:
            return self._geom
        iu = _np.asarray(self.wing_ix)
        jy = _np.asarray(self.wing_jy)
        g = {'ii': iu, 'jj': jy}
        if len(iu) == 0 or len(jy) == 0:
            self._geom = g
            return g
        II, JJ = _np.meshgrid(iu, jy, indexing='ij')
        g['II'], g['JJ'] = II, JJ
        if self.zsurf is not None:
            zs = _np.asarray(self.zsurf[_np.ix_(iu, jy)], dtype=float)
        else:
            zs = _np.zeros((len(iu), len(jy)))
        if self.wing_mask is not None:
            mask = _np.asarray(self.wing_mask[_np.ix_(iu, jy)], dtype=bool)
        else:
            mask = _np.ones(zs.shape, dtype=bool)
        z = self.z
        nz = self.nz
        # face index: z[f] <= zs < z[f+1] ('left' puts zs == z[f] in the
        # cell ABOVE f, i.e. f-1 -- the stock scheme at zs = 0)
        f = _np.searchsorted(z, zs, side='left') - 1
        f = _np.clip(f, 0, nz - 2)
        dzf = _np.maximum(z[f + 1] - z[f], 1e-12)
        th = _np.clip((zs - z[f]) / dzf, 0.0, 1.0)
        # Two different anchors, and they are NOT the same thing:
        #   fcut/th : the z-FACE the sheet/source is applied on. A face can
        #             only bracket a height one-sidedly, so th is the
        #             position of z_s inside the cell [z[f], z[f+1]] and the
        #             two adjacent face stencils are blended (1-th, th);
        #             th = 1 exactly at a node == the stock scheme, which is
        #             what keeps z_s = 0 bit-identical.
        #   ks      : the node the SURFACE QUANTITIES (Cp, w_plane, Kutta
        #             rows) are read around. That one must be the NEAREST
        #             node, not the node above: the stock loads() reads
        #             planes k0-1, k0-2 below the cut, i.e. one plane either
        #             side of the node, and anchoring it at the node ABOVE
        #             makes an epsilon-sized z_s switch the plane pair
        #             (measured: z_s = 1e-4c cost 77% of CL).
        ks = _np.abs(z[None, None, :] - zs[:, :, None]).argmin(axis=-1)
        g.update({'zs': zs, 'mask': mask, 'f': f, 'fcut': f, 'th': th,
                  'ks': ks})
        if self.zsurf is not None and (np.abs(zs).max() > 0):
            gz = _np.gradient(zs, self.xc[iu], self.yc[jy], axis=(0, 1))
            g['zs_x'], g['zs_y'] = gz[0], gz[1]
        else:
            g['zs_x'] = _np.zeros(zs.shape)
            g['zs_y'] = _np.zeros(zs.shape)
        # per-(i,j) section geometry (defaults: the stock single section)
        g['xc_rel'] = _np.repeat(
            (_np.asarray(self.xc[iu]) / self.chord)[:, None], len(jy), axis=1)
        g['alpha'] = _np.full(zs.shape, self.alpha)
        g['slope_t'] = None
        g['slope_c'] = None
        pl = self.planform
        if pl is not None and 'c_j' in pl:
            cj = _np.asarray(pl['c_j'])[jy]
            xlj = _np.asarray(pl['xle_j'])[jy]
            twj = _np.asarray(pl['tw_j'])[jy]
            cc = _np.where(cj > 0, cj, 1.0)
            g['xc_rel'] = ((_np.asarray(self.xc[iu])[:, None] - xlj[None, :])
                           / cc[None, :])
            g['alpha'] = self.alpha + twj[None, :] * _np.ones(
                (len(iu), 1))
            g['alpha'] = _np.broadcast_to(g['alpha'], zs.shape).copy()
            if 'slope_t' in pl:
                g['slope_t'] = pl['slope_t']      # array OR callable
            if 'slope_c' in pl:
                g['slope_c'] = pl['slope_c']
        self._geom = g
        return g

    def _section_slopes(self):
        """Per-(i,j) thickness and camber slopes (local chord frame)."""
        import numpy as _np
        g = self._wing_geom()
        zs = g['zs']
        xr = g['xc_rel']
        if callable(g['slope_t']):
            Y = _np.broadcast_to(_np.asarray(self.yc[g['jj']])[None, :],
                                 xr.shape)
            st = _np.asarray(g['slope_t'](xr, Y), dtype=float)
            sc = (_np.asarray(g['slope_c'](xr, Y), dtype=float)
                  if callable(g['slope_c']) else _np.zeros(xr.shape))
        elif g['slope_t'] is not None:
            st = _np.asarray(g['slope_t'], dtype=float)
            sc = _np.asarray(g['slope_c'], dtype=float)
        elif self.profile is not None:
            st = _np.broadcast_to(
                _np.asarray(self.profile['slope'](xr[:, 0]))[:, None],
                zs.shape)
            sc = _np.broadcast_to(
                _np.asarray(self.profile['camber'](xr[:, 0]))[:, None],
                zs.shape)
        else:
            st = _np.broadcast_to(
                thickness_slope(xr[:, 0], self.thick)[:, None], zs.shape)
            sc = _np.zeros(zs.shape)
        return st, sc

    def _bc_terms(self):
        """Per-(i,j) transpiration velocities wu, wl (upper/lower)."""
        import numpy as _np
        g = self._wing_geom()
        st, sc = self._section_slopes()
        # LOCAL NORMAL PROJECTION. Tangency on z = h(x,y) is
        #     w = h_x*(1 + u/U) + h_y*v - alpha_local      (V_z = sin a ~ a)
        # The u/U factor is dropped exactly as the stock flat BC drops it
        # (it is O(slope*alpha), second order). The h_y*v term is KEPT: for
        # a dihedral/height-varying surface it is the ONLY route by which
        # the tilt of the surface reaches the lift (it is the vertical-
        # velocity form of the normal projection -- the same statement as
        # V.n = 0, since n ~ (h_x, h_y, -1)). It is O(Gamma^2) in the
        # dihedral angle, which is why dihedral moves CL only at second
        # order (cos Gamma - 1 = -Gamma^2/2).
        wu = st + sc + g['zs_x'] - g['alpha']
        wl = -st + sc + g['zs_x'] - g['alpha']
        # MEASURED AND REJECTED (kept behind WING3D_NP_TANGENT=1):
        # the tangential h_y*v term of the normal projection is formally
        # required, but v is a first-difference estimator across the branch
        # cut (the same estimator A10/A11 measured to be jump-contaminated)
        # and its magnitude on this grid is |h_y*v| = 0.10-0.23 alpha for a
        # 5 deg dihedral -- one to two orders of magnitude LARGER than the
        # physical effect it should carry (VLM: dCL = -0.31% at 5 deg,
        # i.e. 0.003 alpha). With it on, mode 7 walks to a negative
        # circulation (measured CL_circ -0.13 vs +0.065 flat). It is
        # therefore OFF by default: a term whose discretisation noise is
        # 30x its physical content is not a model, it is a random number
        # generator with a sign convention.
        if (self.zsurf is not None) and bool((g['zs_y'] != 0).any()) \
                and os.environ.get('WING3D_NP_TANGENT', '0') == '1':
            vu, vl = self._v_surface()
            wu = wu + g['zs_y'] * vu
            wl = wl + g['zs_y'] * vl
        return wu, wl

    def _v_surface(self):
        """Spanwise disturbance velocity v = dphi/dy extrapolated to z_s.

        Same construction as `w_plane` (a first-order difference across the
        cut) but in y: v_s = dphi/dy at the two stencil planes each side,
        linearly extrapolated in z to the surface height.
        """
        import numpy as _np
        g = self._wing_geom()
        ks = g['ks']
        zs = g['zs']
        ii, jj = g['ii'], g['jj']
        nz = self.nz
        # v on every z-plane (central in y; one-sided rows stay 0)
        VY = np.zeros((self.nx, self.ny, nz))
        if self.ny > 2:
            VY[:, 1:-1, :] = (self.phi[:, 2:, :] - self.phi[:, :-2, :]) / (
                self.y[2:] - self.y[:-2])[None, :, None]

        def side(k1, k2):
            k1 = np.clip(k1, 0, nz - 1)
            k2 = np.clip(k2, 0, nz - 1)
            z1, z2 = self.z[k1], self.z[k2]
            den = np.where(np.abs(z2 - z1) > 1e-12, z2 - z1, 1e-12)
            w = (z2 - zs) / den
            a = VY[np.ix_(ii, jj) + (k1,)]
            b = VY[np.ix_(ii, jj) + (k2,)]
            return w * a + (1.0 - w) * b

        ku1 = _np.clip(ks + 1, 0, nz - 1)
        ku2 = _np.clip(ks + 2, 0, nz - 1)
        kl1 = _np.clip(ks - 1, 0, nz - 1)
        kl2 = _np.clip(ks - 2, 0, nz - 1)
        return side(ku1, ku2), side(kl1, kl2)

    def _wake_z(self):
        """Height of the trailing-edge line z_TE(y) per span station.

        The wake sheet must leave the wing AT THE TE HEIGHT of each span
        station. Leaving it on the planar face k0-1 while the wing sheet
        steps up with the dihedral breaks Helmholtz (measured: CL lost
        90% at 2 deg dihedral before this was fixed -- the trailed jump no
        longer equalled the bound TE jump).
        """
        g = self._wing_geom()
        zs = g['zs']
        m = g['mask']
        ng = zs.shape[1]
        zte = np.zeros(ng)
        for q in range(ng):
            rows = np.nonzero(m[:, q])[0]
            zte[q] = zs[rows[-1], q] if rows.size else 0.0
        return zte

    def _wake_faces(self):
        """(face, weight) pairs of the wake sheet per (wake_ix, wake_jy).

        Same in-cell blend as the wing sheet: the trailed jump is applied
        on the two faces bracketing z_TE(y), weighted (1-theta, theta).
        """
        g = self._wing_geom()
        jw = np.asarray(self.wake_jy)
        jy = np.asarray(g['jj'])
        zte = self._wake_z()
        zj = np.interp(np.asarray(self.yc[jw], float),
                       np.asarray(self.yc[jy], float), zte)
        ni = len(self.wake_ix)
        Z = np.repeat(zj[None, :], ni, axis=0)
        f = np.clip(np.searchsorted(self.z, Z, side='left') - 1, 0,
                    self.nz - 2)
        dzf = np.maximum(self.z[f + 1] - self.z[f], 1e-12)
        th = np.clip((Z - self.z[f]) / dzf, 0.0, 1.0)
        return f, f - 1, th, 1.0 - th

    def _stencil_weights(self):
        """(face index, weight) pairs of the two stencils per cell.

        Returns g1, g2 (ni,ng) int face indices and w1, w2 their weights,
        with the mask folded in (inactive cells get weight 0 everywhere).
        """
        import numpy as _np
        g = self._wing_geom()
        f = g['f']
        th = g['th']
        act = g['mask'] & (f >= 1) & (f <= self.nz - 2)
        w1 = _np.where(act, th, 0.0)
        w2 = _np.where(act, 1.0 - th, 0.0)
        return f, f - 1, w1, w2

    def _gradx_planes(self):
        """dphi/dx on every z-plane, shape (nx, ny, nz) (central in x)."""
        import numpy as _np
        G = _np.zeros((self.nx, self.ny, self.nz))
        G[1:-1] = (self.phi[2:, :, :] - self.phi[:-2, :, :]) / (
            self.x[2:] - self.x[:-2])[:, None, None]
        return G

    def _surf_cp_ij(self, G=None, cp_from_g=None):
        """Cp on the upper/lower surface at z_s, per (i,j) cell.

        Linear extrapolation in z to the LOCAL surface height, from the two
        clean planes each side of the snapped sheet (ks+1, ks+2 above and
        ks-1, ks-2 below), exactly as `loads` does for the plane z = 0.
        """
        import numpy as _np
        g = self._wing_geom()
        if G is None:
            G = self._gradx_planes()
        ii, jj = g['ii'], g['jj']
        ks, zs = g['ks'], g['zs']
        nz = self.nz

        def side(k1, k2, sign):
            k1 = _np.clip(k1, 0, nz - 1)
            k2 = _np.clip(k2, 0, nz - 1)
            z1, z2 = self.z[k1], self.z[k2]
            den = _np.where(_np.abs(z2 - z1) > 1e-12, z2 - z1, 1e-12)
            w = (z2 - zs) / den
            a = G[_np.ix_(ii, jj) + (k1,)]
            b = G[_np.ix_(ii, jj) + (k2,)]
            return w * a + (1.0 - w) * b

        gu = side(ks + 1, ks + 2, 1)
        gl = side(ks - 1, ks - 2, -1)
        if cp_from_g is None:
            return -2.0 * gu, -2.0 * gl
        return cp_from_g(gu), cp_from_g(gl)

    def _coeff_A(self, phix):
        gm = self.gamma
        return (1 - self.Minf ** 2) - (gm + 1) * self.Minf ** 2 * phix

    def _face_rho(self, phix_f, phiy_f=None, phiz_f=None):
        """Full-potential face density from face perturbation gradients.

        q = (1+phix, phiy, phiz) (freestream (1,0,0)); rho/rho_inf =
        [1 + (g-1)/2 M^2 (1 - |q|^2)]^(1/(g-1)). q^2 clipped at the
        vacuum-limit value (standard density clipping: prevents NaN
        from fractional powers during transients).
        Transverse gradients averaged to x-faces from nodes when given.
        """
        gm = self.gamma
        q2 = (1.0 + phix_f) ** 2
        if phiy_f is not None:
            q2 = q2 + 0.25 * (phiy_f[:-1, :, :] + phiy_f[1:, :, :]) ** 2 \
                if phiy_f.shape[0] == phix_f.shape[0] + 1 else q2
        if phiz_f is not None:
            q2 = q2 + 0.25 * (phiz_f[:-1, :, :] + phiz_f[1:, :, :]) ** 2 \
                if phiz_f.shape[0] == phix_f.shape[0] + 1 else q2
        q2max = 1 - (0.2 ** (gm - 1) - 1.0) * 2 / max(
            (gm - 1) * self.Minf ** 2, 1e-9)
        q2 = np.minimum(np.maximum(q2, 0.0), q2max)
        rho = (1 + 0.5 * (gm - 1) * self.Minf ** 2 * (1 - q2)) ** (
            1 / (gm - 1))
        rho = np.maximum(rho, 0.2)
        # ---- non-isentropic (entropy) correction -------------------
        # Potential flow is irrotational, so it cannot carry the entropy
        # (and hence vorticity) produced by a shock.  Add the dominant
        # rotational effect: downstream of a shock the total pressure has
        # dropped, so the density is lower than the isentropic value.
        # Detection is streamwise (x) since the wing is x-aligned:
        # track the running max Mach; where the local Mach has fallen
        # below it (i.e. the flow has passed a shock) apply the normal
        # shock total-pressure ratio p02/p01 evaluated at that peak Mach.
        if os.environ.get('WING3D_ENTROPY', '0') == '1':
            a2 = np.maximum(rho, 1e-9) ** (gm - 1)
            M2 = q2 * self.Minf ** 2 / a2
            M2max = np.maximum.accumulate(M2, axis=0)
            shocked = (M2max > 1.0) & (M2 < 0.98 * M2max)
            Mm = np.sqrt(np.maximum(M2max, 1.0))
            num = (gm + 1.0) * Mm ** 2
            den = (gm - 1.0) * Mm ** 2 + 2.0
            p02 = (num / den) ** (gm / (gm - 1.0)) * (
                (gm + 1.0) / (2.0 * gm * Mm ** 2 - (gm - 1.0))
            ) ** (1.0 / (gm - 1.0))
            fac = np.where(shocked, p02, 1.0)
            rho = rho * fac ** (1.0 / gm)
        return rho

    def _xflux(self, phi, dx):
        """x-face fluxes + linearization coefficient + supersonic mask.

        Modes: linear (A const), TSD (A(phi)), fullpot (rho*qx).
        Returns (Fx, A_f, sub).
        """
        phix_f = (phi[1:, :, :] - phi[:-1, :, :]) / dx[:, None, None]
        if self.linear:
            A_f = np.full_like(phix_f, max(1 - self.Minf ** 2, 1e-6))
            return A_f * phix_f, A_f, np.zeros_like(A_f, dtype=bool)
        if getattr(self, 'fullpot', False):
            rho_f = self._frozenRho if self._frozenRho is not None \
                else self._face_rho(phix_f)
            Fx = rho_f * (1.0 + phix_f)
            # face Mach^2: |q|^2 * Minf^2 * (a_inf/a)^2 with
            # (a_inf/a)^2 = rho^-(gamma-1) (isentropic, normalized)
            gm = self.gamma
            q2 = (1.0 + phix_f) ** 2
            Mface = (q2 * self.Minf ** 2
                     / np.maximum(rho_f, 1e-9) ** (gm - 1))
            sub = Mface > 1.0
            if self._frozen is not None:
                sub = self._frozen
            Fup = np.zeros_like(Fx)
            rho_up = rho_f.copy()
            # upwind density: backward face value
            rho_up[1:] = rho_f[:-1]
            Fup[1:] = rho_up[1:] * (
                1.0 + (phi[1:-1, :, :] - phi[:-2, :, :]) /
                dx[:-1, None, None])
            Fx = np.where(sub, Fup, Fx)
            return Fx, rho_f, sub
        if self._frozenA is not None:
            A_f = self._frozenA
        else:
            A_f = self._coeff_A(phix_f)
        Fx = A_f * phix_f
        if self._frozen is not None:
            sub = self._frozen
        else:
            sub = self._sup_mask()
        Fup = np.zeros_like(Fx)
        Fup[1:] = A_f[1:] * (phi[1:-1, :, :] - phi[:-2, :, :]) / \
            dx[:-1, None, None]
        Fx = np.where(sub, Fup, Fx)
        return Fx, A_f, sub

    def _face_coeff(self):
        """Linearization coefficient on x-faces (for frozen epochs and
        the line solver): A (TSD/linear) or rho (full potential)."""
        dx = self.dx
        phix_f = (self.phi[1:, :, :] - self.phi[:-1, :, :]) / dx[:, None, None]
        if self.linear:
            return np.full_like(phix_f, max(1 - self.Minf ** 2, 1e-6))
        if getattr(self, 'fullpot', False):
            if self._frozenRho is not None:
                return self._frozenRho
            return self._face_rho(phix_f)
        return self._coeff_A(phix_f)

    def _nodal_rho(self):
        """Isentropic nodal density from the current phi (fullpot).

        Extracted from residual() so a frozen epoch can reuse a lagged
        copy (self._frozenRhon) and keep the epoch operator affine.
        """
        phi = self.phi
        gx = np.zeros_like(phi)
        gx[1:-1] = (phi[2:] - phi[:-2]) / (
            self.x[2:] - self.x[:-2])[:, None, None]
        gy = np.zeros_like(phi)
        gy[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (
            self.y[2:] - self.y[:-2])[None, :, None]
        gz = np.zeros_like(phi)
        gz[:, :, 1:-1] = (phi[:, :, 2:] - phi[:, :, :-2]) / (
            self.z[2:] - self.z[:-2])[None, None, :]
        q2n = (1 + gx) ** 2 + gy ** 2 + gz ** 2
        gm = self.gamma
        q2max = 1 - (0.2 ** (gm - 1) - 1.0) * 2 / max(
            (gm - 1) * self.Minf ** 2, 1e-9)
        q2n = np.minimum(np.maximum(q2n, 0.0), q2max)
        return np.maximum(
            (1 + 0.5 * (gm - 1) * self.Minf ** 2 * (1 - q2n)) ** (
                1 / (gm - 1)), 0.2)

    def _sup_mask(self):
        """Current supersonic-face mask (for frozen-coefficient epochs)."""
        dx = self.dx
        phix_f = (self.phi[1:, :, :] - self.phi[:-1, :, :]) / dx[:, None, None]
        if getattr(self, 'fullpot', False):
            gm = self.gamma
            q2 = (1.0 + phix_f) ** 2
            rho_f = self._face_rho(phix_f)
            Mface = (q2 * self.Minf ** 2
                     / np.maximum(rho_f, 1e-9) ** (gm - 1))
            return Mface > 1.0
        return self._coeff_A(phix_f) < 0

    def residual(self):
        """Conservative residual, compact finite-volume form.

        Modes: linear (A const), TSD (A(phi), Murman-Cole), fullpot
        (density-coupled exact potential with MC upwinding).
        Node divergence straight from face fluxes (no face-averaging:
        averaging widens the stencil to i+-2 with a checkerboard null
        space that blows up SOR).
        """
        phi, nx, ny, nz = self.phi, self.nx, self.ny, self.nz
        dx, dy, dz = self.dx, self.dy, self.dz
        fp = getattr(self, 'fullpot', False) and not self.linear
        R = np.zeros_like(phi)
        # ---- x ----
        Fx, A_f, _sub = self._xflux(phi, dx)
        # control-volume widths (node-centered)
        wx = np.empty(nx)
        wx[1:-1] = 0.5 * (dx[:-1] + dx[1:])
        wx[0] = dx[0]
        wx[-1] = dx[-1]
        # faces 0..nx-2; node i uses faces i-1, i
        R[1:-1, :, :] += (Fx[1:, :, :] - Fx[:-1, :, :]) / \
            wx[1:-1, None, None]
        # nodal density (fullpot only; TSD/linear use unit density)
        rho_n = None
        if fp:
            rho_n = (self._frozenRhon if self._frozenRhon is not None
                     else self._nodal_rho())
        # ---- y (central, density-weighted in fullpot) ----
        phiy_f = (phi[:, 1:, :] - phi[:, :-1, :]) / dy[None, :, None]
        if fp:
            rho_yf = 0.5 * (rho_n[:, :-1, :] + rho_n[:, 1:, :])
            phiy_f = rho_yf * phiy_f
        wy = np.empty(ny)
        wy[1:-1] = 0.5 * (dy[:-1] + dy[1:])
        wy[0] = dy[0]
        wy[-1] = dy[-1]
        R[:, 1:-1, :] += (phiy_f[:, 1:, :] - phiy_f[:, :-1, :]) / \
            wy[None, 1:-1, None]
        # ---- z (central) + wake branch cut ----
        phiz_f = (phi[:, :, 1:] - phi[:, :, :-1]) / dz[None, None, :]
        k0 = self.k0        # span-cell amplitude k (free bound vortex; ones = prescribed):
        # map each wake span cell onto the wing span grid
        _wpos = np.clip(np.searchsorted(self.wing_jy, self.wake_jy),
                        0, max(len(self.wing_jy) - 1, 0))
        _kw = self.k[_wpos] if len(self.wing_jy) else np.ones(1)
        if len(self.wake_ix) and len(self.wake_jy) and k0 >= 1:
            if self.kutta_on:
                # free circulation: the wake sheet trails exactly the TE
                # value of the bound sheet (Helmholtz), constant in x --
                # no farfield (_circff) blending, which would re-pin it.
                jy = np.asarray(self.wing_jy)
                if self.kutta_morino:
                    # MORINO wake condition: the wake doublet IS the
                    # potential jump the field carries at the TE
                    # (Gamma_wake = phi_up(TE) - phi_lo(TE)), not a
                    # prescribed copy of gamma_b. Helmholtz then holds by
                    # construction and the trailed strength is measured
                    # from the solve instead of being dictated.
                    _Jte = self.te_jump()
                    self.wake_G = np.interp(
                        self.yc, self.yc[jy], _Jte[jy], left=0.0, right=0.0)
                elif self.kutta_free_mu:
                    # mode 5 Helmholtz: the wake trails the TE doublet
                    # mu(x_TE, y) -- one value per span station.
                    self.wake_G = np.interp(
                        self.yc, self.yc[jy],
                        self.mu_b[-1, :], left=0.0, right=0.0)
                else:
                    self.wake_G = np.interp(
                        self.yc, self.yc[jy],
                        self.gamma_b, left=0.0, right=0.0)
            WI, WJ = np.meshgrid(self.wake_ix, self.wake_jy, indexing='ij')
            # muw (wake_G) is circulation PER UNIT SPAN gamma(y)
            # (CL = 2*int(gamma)dy/span, see tests3d/diag_camber.py).
            # The branch-cut potential jump is the SPAN-CELL circulation
            # Gamma = gamma*dy[j]; omitting dy makes the jump span-scaling
            # wrong (measured: CL flips sign as span grows at M0.3).
            _dyj = self.dy[WJ] if os.environ.get(
                'WING3D_JUMP_DYW', '0') == '1' else 1.0
            _g = self.wake_scale * _kw[None, :] * self.wake_G[WJ] * _dyj
            if getattr(self, '_circff', None) is not None \
                    and not self.kutta_on:
                # TSFOIL-style: wake jump interpolated TE -> farfield
                # (box consistency as Gamma evolves)
                xw = self.xc[self.wake_ix]
                t = np.clip((xw - self.chord) / max(
                    self.x[-1] - self.chord, 1e-9), 0.0, 1.0)[:, None]
                gff = self._circff[np.clip(
                    WJ, 0, len(self._circff) - 1)]
                _g = (1 - t) * _g + t * self.wake_scale * gff * (
                    self.dy[WJ] if os.environ.get('WING3D_JUMP_DYW',
                                                  '1') == '1' else 1.0)
            if not self.nonplanar:
                phiz_f[WI, WJ, k0 - 1] -= (self.wake_sign * _g
                                           / self.dz[k0 - 1])
            elif os.environ.get('WING3D_NP_WAKE', '1') != '1':
                phiz_f[WI, WJ, k0 - 1] -= (self.wake_sign * _g
                                           / self.dz[k0 - 1])
            else:
                # Trail the wake at the TE height of each span station.
                # Leaving it on the planar face k0-1 while the wing sheet
                # steps up with the dihedral breaks Helmholtz and the
                # bound circulation collapses (measured, see _wake_z).
                f1, f2, wf1, wf2 = self._wake_faces()
                _wg = self.wake_sign * _g
                for gg, wgt in ((f1, wf1), (f2, wf2)):
                    act = wgt != 0.0
                    if not act.any():
                        continue
                    for gf in np.unique(gg[act]):
                        m = act & (gg == gf)
                        kk = int(gf)
                        phiz_f[WI[m], WJ[m], kk] -= (wgt[m] * _wg[m]
                                                     / self.dz[kk])
        if (self.bound_ramp or self.kutta_bound) and len(self.wing_ix) \
                and len(self.wing_jy) and k0 >= 1:
            BI, BJ = np.meshgrid(self.wing_ix, self.wing_jy, indexing='ij')
            BK = np.tile(np.arange(len(self.wing_jy)),
                         (len(self.wing_ix), 1))
            xc = np.clip(self.xc[self.wing_ix] / self.chord, 0.0, 1.0)
            if self.kutta_free_mu:
                # mode 5: free chordwise doublet distribution. mu_b is
                # the potential jump at the nc control stations; map it
                # onto every wing chordwise node (piecewise linear in x,
                # clamped beyond the end stations).
                jump = np.empty((len(self.wing_ix), len(self.wing_jy)))
                for _p in range(len(self.wing_jy)):
                    jump[:, _p] = np.interp(xc, self.x_ctrl,
                                            self.mu_b[:, _p])
            elif self.kutta_bound:
                # free bound vortex: potential jump = F(x)*gamma_b(y),
                # F(0)=0 at the LE, F(1)=1 at the TE (= the wake jump).
                jump = flat_plate_jump(xc)[:, None] * self.gamma_b[None, :]
            elif self.wing_jump is not None:
                jump = np.asarray(self.wing_jump)
            else:
                ramp = xc[:, None]  # 0 at LE -> 1 at TE
                jump = ramp * self.wake_G[BJ]
            _dyw = self.dy[BJ] if os.environ.get('WING3D_JUMP_DYW',
                                                 '0') == '1' else 1.0
            _js = float(os.environ.get('WING3D_JUMP_SCALE', '1'))
            _src = (self.wake_sign * self.wake_scale * self.k[BK] * jump
                    * _dyw * _js)
            if not self.nonplanar:
                phiz_f[BI, BJ, k0 - 1] -= (_src / self.dz[k0 - 1])
            elif os.environ.get('WING3D_NP_CUT', '1') != '1':
                # ABLATION: keep the cut (potential jump) on the planar
                # face k0-1 even when the body surface is not planar. Only
                # the transpiration source and the surface evaluation move.
                phiz_f[BI, BJ, k0 - 1] -= (_src / self.dz[k0 - 1])
            else:
                # per-(i,j) branch cut on the face that brackets z_s
                g1, g2, w1, w2 = self._stencil_weights()
                for gg, wgt in ((g1, w1), (g2, w2)):
                    act = wgt != 0.0
                    if not act.any():
                        continue
                    for gf in np.unique(gg[act]):
                        m = act & (gg == gf)
                        kk = int(gf)
                        phiz_f[BI[m], BJ[m], kk] -= (wgt[m] * _src[m]
                                                     / self.dz[kk])
        wz = np.empty(nz)
        wz[1:-1] = 0.5 * (dz[:-1] + dz[1:])
        wz[0] = dz[0]
        wz[-1] = dz[-1]
        if fp:
            rho_zf = 0.5 * (rho_n[:, :, :-1] + rho_n[:, :, 1:])
            phiz_f = rho_zf * phiz_f
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
        Upper surface: wu = +dyt/dx + dyc/dx - alpha
        Lower surface: wl = -dyt/dx + dyc/dx - alpha
        (yt = half-thickness, yc = camber/mean line; the camber slope is
        shared by both surfaces, thickness slope is antisymmetric.)
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
        if self.nonplanar:
            # per-(i,j) surface: local slopes, local normal projection
            wu, wl = self._bc_terms()
        else:
            if self.profile is not None:
                slope_t = self.profile['slope'](xc)[:, None]
                slope_c = self.profile['camber'](xc)[:, None]
            else:
                slope_t = thickness_slope(xc, self.thick)[:, None]
                slope_c = np.zeros_like(slope_t)
            # NOTE: no cap here (cap kills LE suction); Newton-Krylov handles
            # the stiff LE source robustly where SOR blew up.
            if self.kutta_thick_only:
                # FREE-CIRCULATION MODE 1 (thickness-only tangency):
                # transpiration is the THICKNESS problem only (+dyt/dx upper,
                # -dyt/dx lower). The lifting problem (camber + alpha) is
                # carried entirely by the bound vortex sheet gamma_b and
                # selected by the Kutta condition. Keeping `slope_c - alpha`
                # here would prescribe part of the lift outside the
                # circulation and make it unreducible.
                wu = slope_t
                wl = -slope_t
            else:
                wu = slope_t + slope_c - self.alpha
                wl = -slope_t + slope_c - self.alpha
        if self.kutta_thick_only and self.nonplanar:
            # same split as the planar branch: the bound sheet carries the
            # lifting part, the source carries the thickness problem only.
            st, _sc = self._section_slopes()
            wu = st
            wl = -st
            g = self._wing_geom()
            if bool((g['zs_y'] != 0).any()) and os.environ.get(
                    'WING3D_NP_TANGENT', '0') == '1':
                vu, vl = self._v_surface()
                wu = wu + g['zs_y'] * vu
                wl = wl + g['zs_y'] * vl
        if self.bl_on and self.bl_dstar_u is not None:
            # Viscous-inviscid coupling: displacement surface
            #     y_u = y_c + y_t + delta*_u ,  y_l = y_c - y_t - delta*_l
            # so its tangency velocity is u_e*d(delta*)/dx with the SAME
            # sign as the thickness slope on each side (a displacement
            # thickness thickens BOTH sides). A14 writes that in the
            # (source/2, vortex) split this BC solves in:
            #     half = (ddu + ddl)/2   -> source/2 : +above, -below
            #                               -> CANCELS in wav (correct: a
            #                                  symmetric thickening produces
            #                                  no mean-plane downwash)
            #     mean = (ddu - ddl)/2   -> viscous CAMBER : +above, +below
            #                               -> the only BL term wav can see
            # A13 read the antisymmetric entry as a bug; the sign is right
            # (the lower surface moves -z), and the non-cancelling part is
            # the asymmetric thickening. Both are exposed by
            # bl_wav_terms() and measured in tests3d/tsd_a14_bl.py.
            ddu, ddl = self._bl_slope()
            _half = 0.5 * (ddu + ddl)      # source/2  (cancels in wav)
            _mean = 0.5 * (ddu - ddl)      # camber    (seen by wav)
            wu = wu + (_half + _mean)
            wl = wl + (-_half + _mean)
        if not self.nonplanar:
            R[II, JJ, k0] += wu / dzp
            R[II, JJ, k0 - 1] += -wl / dzm
            return R
        # non-planar: the sheet sits on the z-face that brackets z_s(x_i,y_j)
        # for THAT (i,j); the two adjacent face stencils are blended by
        # `_stencil_weights` (theta = 1 at a node == the stock scheme).
        g1, g2, w1, w2 = self._stencil_weights()
        dz = self.dz
        for gg, wgt in ((g1, w1), (g2, w2)):
            act = wgt != 0.0
            if not act.any():
                continue
            for gf in np.unique(gg[act]):
                m = act & (gg == gf)
                kk = int(gf)
                R[II[m], JJ[m], kk + 1] += (wgt[m] * wu[m]) / dz[kk + 1]
                R[II[m], JJ[m], kk] += -(wgt[m] * wl[m]) / dz[kk]
        return R

    # ------------------------------------------------------------------
    # viscous-inviscid interaction (integral boundary layer)
    # ------------------------------------------------------------------
    def _surf_edge(self):
        """Edge speed / edge Mach on the wing plane, upper and lower.

        Uses the same linear extrapolation to z=0 as loads()/fp_loads().
        Returns dict with q ('u','l') = u_e/U_inf and me ('u','l') =
        edge Mach, both shaped (len(wing_ix), len(wing_jy)).
        """
        k0 = self.k0

        def gradx(k):
            g = np.zeros((self.nx, self.ny))
            g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
                self.x[2:] - self.x[:-2])[:, None]
            return g

        def surf(k1, k2):
            z1, z2 = self.z[k1], self.z[k2]
            w = abs(z2) / max(abs(z2 - z1), 1e-12)
            return w * gradx(k1) + (1 - w) * gradx(k2)

        if self.cp_decontam():
            # A17: same decontamination as loads()/fp_loads(): bounded
            # upper extrapolation + the exact jump derivative for the
            # lower side, so the BL edge Mach is dz-stable too.
            gxu = self._gx_upper_local()
            gj = self._jump_cell_loading()
            gxl = gxu.copy()
            if gj is not None:
                BI = np.ix_(np.asarray(self.wing_ix),
                            np.asarray(self.wing_jy))
                gxl[BI] = gxu[BI] - gj
        else:
            gxu = surf(min(k0 + 1, self.nz - 1), min(k0 + 2, self.nz - 1))
            gxl = surf(max(k0 - 1, 0), max(k0 - 2, 0))
        iu, jy = np.ix_(np.asarray(self.wing_ix), np.asarray(self.wing_jy))
        gm = self.gamma
        out = {}
        for tag, g in (('u', gxu), ('l', gxl)):
            q = np.maximum(1.0 + g[iu, jy], 1e-3)
            # isentropic: T/T_inf = 1 + (g-1)/2 M^2 (1 - q^2)
            tt = 1.0 + 0.5 * (gm - 1) * self.Minf ** 2 * (1.0 - q ** 2)
            out[tag] = q
            out['me_' + tag] = self.Minf * q / np.sqrt(np.maximum(tt, 0.2))
        return out

    # ------------------------------------------------------------------
    # wing-plane normal velocity (the tangency condition's left side)
    # ------------------------------------------------------------------
    def _gradz_plane(self, k):
        """dphi/dz on z-plane k (central in z; one-sided at the box faces)."""
        g = np.zeros((self.nx, self.ny))
        if 0 < k < self.nz - 1:
            g[:] = (self.phi[:, :, k + 1] - self.phi[:, :, k - 1]) / (
                self.z[k + 1] - self.z[k - 1])
        elif k == 0:
            g[:] = (self.phi[:, :, 1] - self.phi[:, :, 0]) / self.dz[0]
        else:
            g[:] = (self.phi[:, :, -1] - self.phi[:, :, -2]) / self.dz[-1]
        return g

    def w_plane(self, wide=False):
        """Field normal velocity w = dphi/dz at the wing plane, z -> 0.

        MEASURED (tests3d/tsd_a10_probe*.py, g80x32x32, M0.3, AR1.6):
        the branch cut carries a potential jump which is stored as a
        +-J/2 offset on the two nodal planes alongside it. Any 3-point
        stencil that straddles the cut therefore reads jump/dz, not a
        velocity: the wide extrapolation (planes k0+1,k0+2 / k0-1,k0-2)
        gives w = 24.4 * gamma -- 1/dz_core, the plumbing, not the flow.

        The symmetric FIRST-ORDER pair straddling the cut cancels the
        +-J/2 offset to first order and leaves the smooth field: measured
        w = -0.0028 - 0.810*gamma at the 0.75c control point, which is the
        physical downwash scale (thin-airfoil: w ~ -alpha at the chord for
        a lifting distribution). That is the estimator this returns:

            wu = (phi[k0+1] - phi[k0]) / (z[k0+1] - z[k0])
            wl = (phi[k0-1] - phi[k0-2]) / (z[k0-1] - z[k0-2])
            wav = 0.5*(wu + wl)     <- the mean-plane normal velocity

        (wu - wl)/2 is the thickness (antisymmetric) part, which the
        transpiration source already carries. Returns (wu, wl, wav),
        each shaped (nx, ny). `wide=True` gives the old jump-contaminated
        extrapolation, kept for diagnostics only.
        """
        k0 = self.k0
        if wide:
            def surf(k1, k2):
                z1, z2 = self.z[k1], self.z[k2]
                w = abs(z2) / max(abs(z2 - z1), 1e-12)
                return w * self._gradz_plane(k1) + (1 - w) * \
                    self._gradz_plane(k2)

            ku1 = min(k0 + 1, self.nz - 1)
            ku2 = min(k0 + 2, self.nz - 1)
            kl1 = max(k0 - 1, 0)
            kl2 = max(k0 - 2, 0)
            wu = surf(ku1, ku2)
            wl = surf(kl1, kl2)
            return wu, wl, 0.5 * (wu + wl)
        kp = min(k0 + 1, self.nz - 1)
        km = max(k0 - 2, 0)
        wu = (self.phi[:, :, kp] - self.phi[:, :, k0]) / max(
            self.z[kp] - self.z[k0], 1e-12)
        wl = (self.phi[:, :, k0 - 1] - self.phi[:, :, km]) / max(
            self.z[k0 - 1] - self.z[km], 1e-12)
        if self.nonplanar:
            # per-(i,j) stencil: the surface plane moves with (x,y)
            g = self._wing_geom()
            ii, jj, ks = g['ii'], g['jj'], g['ks']
            kpu = np.clip(ks + 1, 0, self.nz - 1)
            kmu = np.clip(ks - 2, 0, self.nz - 1)
            klu = np.clip(ks - 1, 0, self.nz - 1)
            zu = np.maximum(self.z[kpu] - self.z[ks], 1e-12)
            zl = np.maximum(self.z[klu] - self.z[kmu], 1e-12)
            ix = np.ix_(ii, jj)
            wu[ix] = (self.phi[ix + (kpu,)] - self.phi[ix + (ks,)]) / zu
            wl[ix] = (self.phi[ix + (klu,)] - self.phi[ix + (kmu,)]) / zl
        return wu, wl, 0.5 * (wu + wl)

    def te_jump(self):
        """Potential jump actually carried by the field at the TE.

        phi(k0) - phi(k0-1) at the last wing chord station, full ny. This
        is the *measured* jump used by the Morino wake condition
        (Gamma_wake = phi_up(TE) - phi_lo(TE)); compare it with the
        prescribed jump F(1)*gamma_b = gamma_b to see whether that
        condition carries any information (measured, see
        tsd_a10_probe.py).
        """
        k0 = self.k0
        ite = self.wing_ix[-1]
        if self.nonplanar:
            g = self._wing_geom()
            p = int(np.searchsorted(np.asarray(self.wing_ix), ite))
            p = min(p, len(self.wing_ix) - 1)
            ks = g['ks'][p]
            jy = np.asarray(self.wing_jy)
            out = np.zeros(self.ny)
            out[jy] = (self.phi[ite, jy, ks] - self.phi[ite, jy, ks - 1])
            return out
        return self.phi[ite, :, k0] - self.phi[ite, :, k0 - 1]

    def bl_march(self, Re=None, xtr=None):
        """March an integral BL along every span strip (both surfaces).

        Edge velocity comes from the inviscid surface gradient, densified
        onto a fine chordwise grid (the TSD wing plane has only ~13
        chordwise stations, too coarse to differentiate Ue directly).
        Laminar Thwaites + Michel/fixed transition + Head entrainment
        turbulent closure (wing3d/ibl.py). Returns delta*/c arrays.
        """
        from wing3d.ibl import solve_streamline
        if Re is None:
            Re = self.bl_Re
        if xtr is None:
            xtr = self.bl_xtr
        ed = self._surf_edge()
        xc = self.xc[np.asarray(self.wing_ix)] / self.chord
        nf = 161
        xf = np.linspace(xc[0], xc[-1], nf)
        nu = 1.0 / max(Re, 1e-30)  # chord units, U_inf = 1
        du = np.zeros((len(xc), len(self.wing_jy)))
        dl = np.zeros_like(du)
        seps = []
        for tag, dst in (('u', du), ('l', dl)):
            q = ed[tag]
            me = ed['me_' + tag]
            for jj in range(q.shape[1]):
                ucol = np.interp(xf, xc, q[:, jj])
                mcol = np.interp(xf, xc, me[:, jj])
                # light smoothing: raw 13-point inviscid slopes are noisy
                # enough to make dUe/dx (Thwaites lambda) oscillate
                k = np.ones(5) / 5.0
                ucol = np.convolve(np.pad(ucol, 2, mode='edge'), k,
                                   mode='valid')
                r = solve_streamline(xf, ucol, nu, x_tr_forced=xtr,
                                     me=mcol)
                dstar = np.clip(r['dstar'], 0.0, 0.05)
                dst[:, jj] = np.interp(xc, xf, dstar)
                if r['separated']:
                    i_sep = int(r['i_sep']) if r['i_sep'] >= 0 else nf - 1
                    seps.append({'side': tag, 'j': int(jj),
                                 'y': float(self.yc[self.wing_jy[jj]]),
                                 'x_sep': float(xf[min(i_sep, nf - 1)]),
                                 'H': float(r['H'][min(i_sep, nf - 1)]),
                                 'bubble': bool(r['bubble'])})
        self.bl_sep = seps
        return {'dstar_u': du, 'dstar_l': dl, 'sep': seps,
                'Re': Re, 'xtr': xtr}

    def _bl_slope(self):
        """u_e * d(delta*)/dx on the wing plane (upper, lower)."""
        xc = self.xc[np.asarray(self.wing_ix)] / self.chord
        ed = self._surf_edge()
        ddu = np.gradient(self.bl_dstar_u, xc, axis=0) * ed['u']
        ddl = np.gradient(self.bl_dstar_l, xc, axis=0) * ed['l']
        # diagnostic gain (like WING3D_JUMP_SCALE): lets the sensitivity
        # dCL/d(displacement) be measured and the coupling checked for
        # liveness. Default 1 = the physical term.
        sc = float(os.environ.get('WING3D_BL_SCALE', '1'))
        return ddu * sc, ddl * sc

    def bl_wav_terms(self):
        """Decompose the BL transpiration into what wav sees and what not.

        Returns a dict with
          ddu, ddl : the two surfaces' u_e*d(delta*)/dx (z-velocity units)
          half     : (ddu+ddl)/2 -- the source/2 part, enters wu=+half,
                     wl=-half, so it CANCELS in wav=(wu+wl)/2
          mean     : (ddu-ddl)/2 -- the viscous camber, enters both sides
                     with +, so it IS the BL's contribution to the tangency
                     row, and hence the only route from viscosity to CL
          wav_alpha: mean/alpha -- its size in units of the angle of attack
        """
        if self.bl_dstar_u is None:
            return {}
        ddu, ddl = self._bl_slope()
        return {'ddu': ddu, 'ddl': ddl,
                'half': 0.5 * (ddu + ddl), 'mean': 0.5 * (ddu - ddl),
                'wav_alpha': 0.5 * (ddu - ddl) / max(abs(self.alpha), 1e-30)}

    def bl_update(self, relax=None, verbose=False):
        """One viscous-inviscid coupling iteration.

        Cp -> integral BL -> relaxed transpiration. Under-relaxed because
        the displacement term is a fixed RHS source (the operator does
        not see it), so the coupling is a plain fixed point and delta*
        is a very stiff unknown on this coarse chordwise grid.
        """
        if relax is None:
            relax = self.bl_relax
        r = self.bl_march()
        dnu, dnl = r['dstar_u'], r['dstar_l']
        if self.bl_dstar_u is None:
            ch = max(float(np.max(np.abs(dnu))), float(np.max(np.abs(dnl))))
            self.bl_dstar_u = dnu.copy()
            self.bl_dstar_l = dnl.copy()
        else:
            ch = max(float(np.max(np.abs(dnu - self.bl_dstar_u))),
                     float(np.max(np.abs(dnl - self.bl_dstar_l))))
            self.bl_dstar_u = (1 - relax) * self.bl_dstar_u + relax * dnu
            self.bl_dstar_l = (1 - relax) * self.bl_dstar_l + relax * dnl
        info = {'dstar_max': float(max(self.bl_dstar_u.max(),
                                       self.bl_dstar_l.max())),
                'change': ch, 'n_sep': len(r['sep']), 'sep': r['sep'],
                'Re': r['Re'], 'xtr': r['xtr']}
        self.bl_hist.append(info)
        if verbose:
            print('    BL: dstar*max/c=%.4e  d=%.2e  sep=%d' % (
                info['dstar_max'], ch, info['n_sep']), flush=True)
        return info

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

    def solve(self, itmax=600, tol=1e-5, verbose=True, ramp=100,
              freeze_every=0, init_linear=True, use_lines=False):
        """Relaxation solver (red-black SOR or x-line relaxation).

        ramp: iterations over which wake_scale goes 0 -> 1 (avoids
        startup shock from impulsive circulation).
        freeze_every: nonlinear mode: refresh the frozen supersonic map
        every N iters (0 = live map). init_linear: solve the linear
        problem first and use it as the nonlinear init (continuation).
        use_lines: x-line relaxation (Thomas per x-line, alternating
        sweep direction) instead of point SOR -- far better on
        stretched grids and required for stiff transonic cases.
        """
        if init_linear and not self.linear:
            # linear init WITH the same wake/farfield (lifting from step 0)
            lin = TSDSolver(self.Minf, alpha_deg=np.degrees(self.alpha),
                            chord=self.chord, span=self.span,
                            thick=self.thick, gamma=self.gamma,
                            nx=self.nx, ny=self.ny, nz=self.nz,
                            omega=1.2,
                            wake_gamma=None, farfield=self.farfield,
                            linear=True, bound_ramp=self.bound_ramp,
                            profile=self.profile)
            lin.wake_G = self.wake_G.copy()
            lin.wake_scale = 1.0
            lin.solve(itmax=400, tol=1e-6, verbose=False, ramp=0,
                      init_linear=False)
            # transplant interior (keep our farfield box values)
            self.phi[1:-1, 1:-1, 1:-1] = lin.phi[1:-1, 1:-1, 1:-1]
            if verbose:
                print('  linear init done', flush=True)
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
        beta2 = max(1 - self.Minf ** 2, 1e-6)
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
            # frozen-coefficient epochs (nonlinear): refresh supersonic
            # map + full A every N iters; SOR converges the fixed
            # linear operator inside each epoch
            if freeze_every and not self.linear:
                if it % freeze_every == 0:
                    self._frozen = self._sup_mask()
                    self._frozenA = self._face_coeff()
                    if getattr(self, 'fullpot', False):
                        dx = self.dx
                        phix = (self.phi[1:, :, :] - self.phi[:-1, :, :]) / \
                            dx[:, None, None]
                        self._frozenRho = self._face_rho(phix)
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
            if use_lines:
                self.line_relax(omega=0.8,
                                direction=+1 if it % 2 == 0 else -1)
                if not np.all(np.isfinite(self.phi)):
                    bad = np.argwhere(~np.isfinite(self.phi))[0]
                    print('NAN at it %d loc %s x=%.3f y=%.3f z=%.3f' % (
                        it, bad, self.x[min(bad[0], self.nx - 1)],
                        self.y[min(bad[1], self.ny - 1)],
                        self.z[min(bad[2], self.nz - 1)]), flush=True)
                    break
            else:
                for mask in (red, ~red):
                    R = self.residual()
                    upd = np.zeros_like(self.phi)
                    # SOR for Ax=b with A_ii<0 (discrete Laplacian-like):
                    # x += (w/A_ii)(b-Ax) = +w*R/D, D = |A_ii| > 0.
                    upd[mask] = +self.omega * R[mask] / D[mask]
                    self._apply_bc(upd)
                    self.phi[mask] += upd[mask]
            if it % 50 == 0 or it == itmax - 1:
                R = self.residual()
                r = np.abs(R[1:-1, 1:-1, 1:-1]).max()
                if verbose:
                    nsup = int(self._sup_mask().sum()) if not self.linear \
                        else 0
                    print(f'  it {it}: max|R|={r:.3e} nsup={nsup}')
                if r < tol:
                    break
        self._frozen = None  # release map (live evaluation afterwards)
        self._frozenA = None
        self._frozenRho = None
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

    # ------------------------------------------------------------------
    # branch-cut-decontaminated surface estimator (A17)
    # ------------------------------------------------------------------
    def cp_decontam(self):
        """True when the decontaminated surface estimator is enabled.

        WING3D_CPDECONTAM=1. Default OFF reproduces the stock estimator
        bit-identically (no float operation is inserted anywhere).

        WHY (measured, tests3d/a17_cpdecontam.py, g40x24x16 M0.3 AR1.6):
        the stock Cp comes from `gradx` differenced on the planes
        k0+-1, k0+-2 and linearly extrapolated to z=0 with the weight
        w = |z2|/|z2-z1| > 1. That weight AMPLIFIES the near plane by
        1 + dz1/dz2, so the reported Cp jump scales with the z
        distribution: CL = 0.0703/0.0554/0.0536/0.0438 for nz =
        12/16/24/32 (47.5% spread) and a rigid 0.25/0.50 dz offset moves
        it -9.7%/+25.0%, while the circulation the mu sheet carries
        (2*int Gamma dy/S) is 0.05156 at every nz (0.0% spread). The
        whole spread is therefore in the ESTIMATOR, not the solution.

        What is NOT true: the stored phi does not carry a +-J/2 offset.
        Measured, phi[k0] - phi[k0-1] equals the prescribed jump J to 2%
        at every station and every nz -- the residual already removes the
        J/dz it imposes on the cut face, so the two nodal planes ARE the
        true one-sided surface values. Subtracting a further J/2 (the
        A16 prescription) is not a decontamination: it re-injects the
        whole jump and raises the spread to 62%.

        What IS exact and usable is the sheet identity

            dphi/dx(z=0+) - dphi/dx(z=0-) = dJ/dx        (dCp = 2 dJ/dx)

        with J the jump the field carries (known exactly). The upper
        side is taken from the field with a BOUNDED local extrapolation
        (`_gx_upper_local`, weight <= 1, no dz amplification); the lower
        side is the upper one minus the exact jump derivative
        (`_jump_cell_loading`, cell-integrated so it telescopes to
        Gamma). That removes the z-extrapolation from the lift entirely:
        measured dz spread 47.5% -> 1.1%, z-offset -9.7/+25.0% ->
        -0.1/-0.4%, CL 0.0483-0.0486 vs the sheet's own 0.05156 (-6%).

        Consequence to report, not hide: on the stock 6-node chord the
        stock extrapolation was also (accidentally) compensating the
        unresolved leading-edge suction cell, so the decontaminated CL
        sits 13% BELOW the stock number at nz=16 and ~21% below
        thin-airfoil finite-AR truth: the residual error is CHORDWISE
        resolution (nx), not dz. See tests3d/a17_cpdecontam.py.
        """
        return os.environ.get('WING3D_CPDECONTAM', '0') == '1'

    def _gradx_plane(self, k):
        g = np.zeros((self.nx, self.ny))
        g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
            self.x[2:] - self.x[:-2])[:, None]
        return g

    def _gx_upper_local(self):
        """dphi/dx on the upper side of the sheet, extrapolated to z=0.

        Bounded: node k0 is the upper side of the cut face k0-1 and |z[k0]|
        is at most one cell, so the 2-point weight is |t| <= 1 and no
        dz-dependent amplification enters (contrast the stock weight
        w = |z2|/|z2-z1| > 1 on the two FAR planes).
        """
        k0 = self.k0
        k1 = min(k0 + 1, self.nz - 1)
        g0 = self._gradx_plane(k0)
        if k1 == k0:
            return g0
        g1 = self._gradx_plane(k1)
        dz = self.z[k1] - self.z[k0]
        t = 0.0 if abs(dz) < 1e-14 else -float(self.z[k0]) / dz
        t = min(max(t, -1.0), 1.0)
        return g0 + t * (g1 - g0)

    def sheet_jump(self):
        """Potential jump the STORED field carries on the wing nodes.

        phi[k0] - phi[k0-1] at the wing footprint, shape
        (len(wing_ix), len(wing_jy)). Measured to equal the prescribed
        doublet strength to 2% at every station and every nz -- this is
        the exact quantity the decontamination uses, so it is read from
        the field rather than re-derived from the closure mode.
        """
        k0 = self.k0
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        if k0 < 1 or iu.size == 0 or jy.size == 0:
            return np.zeros((iu.size, jy.size))
        return (self.phi[iu, :, k0] - self.phi[iu, :, k0 - 1])[:, jy]

    def _jump_cell_loading(self):
        """d(jump)/dx per wing chord cell, cell-integrated.

        Faces are the averages of the nodal jump, with the two END faces
        pinned to the physical sheet ends (no jump ahead of the LE, the
        TE jump behind the TE). Then sum_i gamma_i dx_i = J_TE = Gamma by
        construction, so CL = 2 int gamma dx ties to the model's own
        circulation. A plain nodal finite difference loses the LE cell,
        which on the stock 6-node chord carries 44% of Gamma (measured:
        CL 0.0333 instead of 0.0483).
        """
        J = self.sheet_jump()
        n = J.shape[0]
        if n == 0:
            return None
        if n == 1:
            return np.zeros_like(J)
        Jf = np.empty((n + 1, J.shape[1]))
        Jf[0] = 0.0
        Jf[1:n] = 0.5 * (J[:-1] + J[1:])
        Jf[n] = J[-1]
        iu = np.asarray(self.wing_ix)
        dxi = np.empty(n)
        dxi[:-1] = np.diff(self.x[iu])
        dxi[-1] = dxi[-2]
        return (Jf[1:] - Jf[:-1]) / dxi[:, None]

    def _apply_cpdecontam(self, cpu, cpl, cp_from_g=None):
        """Replace the wing-footprint Cp by the decontaminated pair.

        cp_from_g : callable gx -> Cp (None = TSD order, Cp = -2 gx).
        Outside the wing footprint both fields keep the stock values, so
        every other consumer sees an unchanged array.
        """
        k0 = self.k0
        gj = self._jump_cell_loading()
        if k0 < 1 or gj is None:
            return cpu, cpl
        gxu = self._gx_upper_local()
        gxl = gxu.copy()
        BI = np.ix_(np.asarray(self.wing_ix), np.asarray(self.wing_jy))
        gxl[BI] = gxu[BI] - gj
        if cp_from_g is None:
            return -2 * gxu, -2 * gxl
        return cp_from_g(gxu), cp_from_g(gxl)

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
        if self.nonplanar:
            return self._loads_ij()
        cpu = surf_cp(ku1, ku2)
        cpl = surf_cp(kl1, kl2)
        if self.cp_decontam():
            # A17: exact-jump split at the sheet (see _apply_cpdecontam).
            # The non-planar path (_loads_ij) is deliberately NOT changed:
            # it extrapolates per (i,j) to its own surface height.
            cpu, cpl = self._apply_cpdecontam(cpu, cpl)
        CL = CD = 0.0
        for i in self.wing_ix:
            xc = self.xc[i] / self.chord
            if self.profile is not None:
                st = self.profile['slope'](xc)
            else:
                st = thickness_slope(xc, self.thick)
            for j in self.wing_jy:
                dA = dxi[min(i, len(dxi) - 1)] * self.dy[min(j, len(self.dy) - 1)]
                CL += (cpl[i, j] - cpu[i, j]) * dA
                CD += (cpu[i, j] - cpl[i, j]) * st * dA
        S = self.chord * self.span
        return {'CL': CL / S, 'CD': CD / S, 'cpu': cpu, 'cpl': cpl}

    def fp_loads(self):
        """Full-potential loads: isentropic Cp from surface Mach.

        Cp = (2/gM^2)(p/pinf - 1) with p/pinf from local isentropic
        Mach (extrapolated like loads()). Use instead of loads() in
        fullpot mode (Cp = -2 phix is TSD-order, wrong at shocks).
        """
        k0 = self.k0
        dxi = self.dx
        gm = self.gamma

        def gradx(k):
            g = np.zeros((self.nx, self.ny))
            g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
                self.x[2:] - self.x[:-2])[:, None]
            return g

        def surf_gx(k1, k2):
            z1 = self.z[k1]
            z2 = self.z[k2]
            w = abs(z2) / max(abs(z2 - z1), 1e-12)
            return w * gradx(k1) + (1 - w) * gradx(k2)

        ku1 = min(k0 + 1, self.nz - 1)
        ku2 = min(k0 + 2, self.nz - 1)
        kl1 = max(k0 - 1, 0)
        kl2 = max(k0 - 2, 0)
        gxu = surf_gx(ku1, ku2)
        gxl = surf_gx(kl1, kl2)

        def cp_of_g(gx):
            q2 = np.maximum((1.0 + gx) ** 2, 1e-9)
            # local Mach^2 from q (isentropic, normalized)
            # p/pinf via Mach: solve M from q requires total; instead
            # use density-based: rho from q, p/pinf = rho^gm
            rho = np.maximum(
                np.maximum(1 + 0.5 * (gm - 1) * self.Minf ** 2 * (1 - q2),
                           1e-9) ** (1 / (gm - 1)), 0.2)
            p_rat = rho ** gm
            return 2.0 / (gm * self.Minf ** 2) * (p_rat - 1.0)

        if self.nonplanar:
            return self._loads_ij(cp_of_g=cp_of_g)
        cpu = cp_of_g(gxu)
        cpl = cp_of_g(gxl)
        if self.cp_decontam():
            # A17: exact-jump split at the sheet (see _apply_cpdecontam).
            cpu, cpl = self._apply_cpdecontam(cpu, cpl, cp_from_g=cp_of_g)
        CL = CD = 0.0
        for i in self.wing_ix:
            xc = self.xc[i] / self.chord
            if self.profile is not None:
                st = self.profile['slope'](xc)
            else:
                st = thickness_slope(xc, self.thick)
            for j in self.wing_jy:
                dA = dxi[min(i, len(dxi) - 1)] * self.dy[min(j, len(self.dy) - 1)]
                CL += (cpl[i, j] - cpu[i, j]) * dA
                CD += (cpu[i, j] - cpl[i, j]) * st * dA
        S = self.chord * self.span
        return {'CL': CL / S, 'CD': CD / S, 'cpu': cpu, 'cpl': cpl}

    def _cell_mask(self):
        """(nc, ng) mask of mu cells whose planform cell carries a body.

        Cells outside the planform have no physics: their tangency row is
        replaced by `mu = 0` (an identity row), otherwise the system picks
        up rows that are identically zero and the Newton solve is singular.
        The stock rectangle has every cell active, so this is a no-op.
        """
        ng = len(self.wing_jy)
        nc = self.n_ctrl
        if not self.nonplanar or self.wing_mask is None:
            return np.ones((nc, ng), dtype=bool)
        g = self._wing_geom()
        m = np.asarray(g['mask'], bool)
        ii = np.asarray(g['ii'])
        rows = np.array([int(np.argmin(np.abs(ii - iu)))
                         for iu in np.asarray(self.wing_ix)[
                             np.clip(np.asarray(self.i_ctrl, dtype=int), 0,
                                     len(self.wing_ix) - 1)]])
        return m[rows]

    def _zero_rows(self, w, g_):
        """Replace the rows of inactive mu cells by the identity row."""
        cm = self._cell_mask()
        if cm.all():
            return w, g_
        w = np.array(w, dtype=float, copy=True)
        g_ = np.array(g_, dtype=float, copy=True)
        w[~cm] = self.mu_b[~cm]
        g_[~cm] = 0.0
        return w, g_

    def _camber_at(self, xr):
        """Camber slope d(yc)/dx at (control station, span station).

        `xr` is (nc, ng) in LOCAL chord units, so a tapered/twisted
        multi-section planform gets its own camber at each station instead
        of the root section's.
        """
        xr = np.asarray(xr, float)
        pl = self.planform
        if pl is not None and callable(pl.get('slope_c')):
            Y = np.broadcast_to(
                np.asarray(self.yc[np.asarray(self.wing_jy)])[None, :],
                xr.shape)
            return np.asarray(pl['slope_c'](xr, Y), float)
        if self.profile is not None:
            return np.asarray(self.profile['camber'](xr), float)
        return np.zeros(xr.shape)

    def _geom_target(self, slope_c, xcp_abs=None):
        """(nc, ng) tangency target: the mean-plane normal velocity a
        non-planar surface REQUIRES, i.e. dyc/dx + dz_s/dx - alpha_local
        (+ the tangential h_y*v term of the local normal projection).

        Reduces exactly to the stock (slope_c - alpha) broadcast when the
        surface is planar, so nothing changes by default.
        """
        nc = len(slope_c)
        base = np.asarray(slope_c, dtype=float) - self.alpha
        if not self.nonplanar:
            if base.ndim == 1:
                return base[:, None] * np.ones((1, len(self.wing_jy)))
            return base
        g = self._wing_geom()
        ii = np.asarray(g['ii'])
        ng = len(g['jj'])
        if xcp_abs is None:
            rows = np.asarray(self.i_ctrl)
        else:
            rows = np.array([int(np.argmin(np.abs(self.xc[ii] - xq)))
                             for xq in np.asarray(xcp_abs, float)])
        tgt = base + g['zs_x'][rows] - (g['alpha'][rows] - self.alpha)
        if bool((g['zs_y'] != 0).any()) and os.environ.get(
                'WING3D_NP_TANGENT', '0') == '1':
            vu, _vl = self._v_surface()
            tgt = tgt + g['zs_y'][rows] * vu[rows]
        return tgt * np.ones((1, ng))

    def _loads_ij(self, cp_of_g=None):
        """Lift/drag for a general (masked, non-planar) surface.

        Same Cp definition and the same extrapolation-to-the-surface rule
        as `loads`, applied per (i,j) at the local surface height. The area
        weights are the cell areas of the ACTIVE mask, and the reference
        area is the planform area (plan S) when one is installed, else
        chord*span -- so a tapered flying wing is normalised by its own
        area, not by a rectangle.
        """
        cpu_ij, cpl_ij = self._surf_cp_ij(cp_from_g=cp_of_g)
        g = self._wing_geom()
        ii, jj = g['ii'], g['jj']
        st, _sc = self._section_slopes()
        dxi = self.dx
        # NOTE: sequential accumulation in the SAME (i-major) order as the
        # planar `loads` loop, so that zsurf = 0 reproduces CL bit for bit
        # (np.sum's pairwise order does not).
        CL = CD = 0.0
        for p, i in enumerate(ii):
            dAi = dxi[min(int(i), len(dxi) - 1)]
            for q, j in enumerate(jj):
                if not g['mask'][p, q]:
                    continue
                dA = dAi * self.dy[min(int(j), len(self.dy) - 1)]
                CL += (cpl_ij[p, q] - cpu_ij[p, q]) * dA
                CD += (cpu_ij[p, q] - cpl_ij[p, q]) * st[p, q] * dA
        S = self.S_ref()
        out = {'CL': CL / S, 'CD': CD / S}
        # keep the (nx,ny) shape the planar path returns for consumers
        full_u = np.zeros((self.nx, self.ny))
        full_l = np.zeros((self.nx, self.ny))
        full_u[np.ix_(ii, jj)] = cpu_ij
        full_l[np.ix_(ii, jj)] = cpl_ij
        out['cpu'], out['cpl'] = full_u, full_l
        return out

    def S_ref(self):
        """Reference area: the installed planform area, else chord*span."""
        pl = self.planform
        if pl is not None and 'chord' in pl and 'y' in pl:
            y = np.asarray(pl['y'], dtype=float)
            c = np.asarray(pl['chord'], dtype=float)
            S = float(np.trapezoid(c, y))
            if S > 0:
                return S
        return float(self.chord * self.span)

    def kutta_residual(self, upstream=0):
        """TE pressure-jump per wing span cell (Kutta rows).

        Returns vector length len(wing_jy): dCp = cpl - cpu at the
        last wing chord cell (upstream=0) or `upstream` cells ahead
        of TE (upstream>0, off the cut-junction singularity).
        Same extrapolation as loads(). Zero at a Kutta-satisfying
        state; drives the free amplitudes k.
        """
        k0 = self.k0

        def gradx(k):
            g = np.zeros((self.nx, self.ny))
            g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
                self.x[2:] - self.x[:-2])[:, None]
            return g

        def surf_cp(k1, k2):
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
        ite = self.wing_ix[-1 - upstream] if upstream < len(
            self.wing_ix) else self.wing_ix[0]
        return np.array([cpl[ite, j] - cpu[ite, j]
                         for j in self.wing_jy])

    def tangency_residual(self, xc_cp=None):
        """Tangency row: the field's normal velocity at the mean plane.

        For each span cell j the physical condition is that the normal
        velocity at the control station balances the geometric slope:

            wav(x_cp, y_j) - (dyc/dx(x_cp) - alpha) = 0

        wav = 0.5*(wu + wl) is the mean-plane normal velocity
        (``w_plane``: dphi/dz extrapolated to z=0 from the planes k0+-1,
        k0+-2, i.e. read from phi *values*, never from the branch-cut
        face flux).  Its thickness half (wu - wl)/2 is the antisymmetric
        part, which the transpiration source already carries; only the
        symmetric (lifting) part is left for this row, which is why the
        row must be paired with thickness-only transpiration
        (kutta_thick_only): if the transpiration also imposed
        (dyc/dx - alpha) the row would be satisfied by construction and
        could not select gamma.

        Returned length len(wing_jy); scaled by kutta_w in the coupled
        residual. Unlike kutta_residual (which reads phi gradients at the
        TE and therefore only constrains the *local* sheet amplitude),
        wav at station j contains the downwash of the whole trailed
        system, so the row couples all span stations -> elliptic loading.
        """
        k0 = self.k0
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        if len(iu) == 0 or len(jy) == 0:
            return np.zeros(0)
        xc = self.xc[iu] / self.chord
        xcp = float(self.tang_xc if xc_cp is None else xc_cp)
        wu, wl, wav = self.w_plane()
        if self.tang_mode == 'gal':
            # LOADING-WEIGHTED (Galerkin) form of the same condition: the
            # chordwise tangency is reduced with the only shape the sheet
            # has, the flat-plate loading A(x). One point per span cell
            # cannot represent the chordwise problem, and a single 0.75c
            # station reads the LOCAL camber slope there (dfvlr: +2.5% at
            # 0.75c -> zero circulation, measured in a10 ladder A). The
            # weighted form is the standard reduction of the chordwise
            # condition to a section equation.
            A = (2.0 / np.pi) * np.sqrt(
                np.maximum(1.0 - xc, 0.0) / np.maximum(xc, 1e-9))
            wgt = A / max(float(np.trapezoid(A, xc)), 1e-30)
            if self.profile is not None:
                slope_c = np.asarray(self.profile['camber'](xc), dtype=float)
            else:
                slope_c = np.zeros(len(xc))
            g = slope_c - self.alpha
            r = np.zeros(len(jy))
            for p, j in enumerate(jy):
                r[p] = float(np.trapezoid(wgt * (wav[iu, j] - g), xc))
            return r
        icp = int(iu[int(np.argmin(np.abs(xc - xcp)))])
        x_at = self.xc[icp] / self.chord
        if self.profile is not None:
            slope_c = float(self.profile['camber'](x_at))
        else:
            slope_c = 0.0
        w_geom = slope_c - self.alpha
        return np.array([wav[icp, j] for j in jy]) - w_geom

    # ------------------------------------------------------------------
    # mode 6: ANALYTIC DOUBLET-LATTICE INFLUENCE KERNEL
    # ------------------------------------------------------------------
    def _dlm_panels(self):
        """Panel lattice for the mode-6 kernel: wing sheet + trailed wake.

        The mode-5 tangency rows read the normal velocity out of the TSD
        field (`w_plane`), which is a FIRST-ORDER finite difference across
        the branch cut on a grid whose dz_core is 1.5% of the chord. A11
        measured that estimator's own chordwise kernel to be 3.7x
        stronger aft than at the LE for a physically uniform downwash
        (0.45*alpha at the LE cell to 1.7*alpha aft), which flattened mu
        and cost ~35% of the circulation.

        Here the downwash of the doublet sheet is computed ANALYTICALLY
        instead. A constant doublet panel of strength mu is equivalent to
        a vortex ring of circulation mu around its edges, and the velocity
        a vortex ring induces at a point IN ITS OWN PLANE is purely normal
        to that plane -- so the classical doublet-lattice influence
        coefficient is exactly the z-component of the ring kernel. That
        kernel already exists in this codebase as
        `wing3d.ckernels.loop_vel_block` (used by the 3D Morino panel
        solver), so this is a REUSE of the same Biot-Savart source, not a
        new singular integral.

        Geometry: chordwise cells between control STATIONS (cell i spans
        [x_{i-1}, x_i], LE edge = 0), so the cell's unknown loading
        Gamma_i = mu_i - mu_{i-1} sits on the chord it actually covers;
        span cells are centred on a station (edges at the midpoints, tips
        at +-span/2). One horseshoe per cell: bound filament at the 1/4
        chord, legs trailing downstream, control point at the 3/4 chord
        (classical 1/4-3/4). The trailed wake is therefore inside the
        kernel, so the downwash it produces is part of the tangency
        equation.

        MEASURED alternatives (both rejected, numbers in
        tests3d/a12_kutta_report.md): a constant-doublet panel reaching
        x=0 puts the whole mu(0.043c) jump on the LE edge and reads
        +0.41*alpha at the LE row (the exact solution needs -alpha); a
        sub-panel stack (SUBX>1) does not converge, because the doublet
        kernel is 1/r^2 on the sheet and a midpoint stack of constant
        sub-panels gives +22*alpha mid-chord instead of -alpha.
        """
        ys = np.asarray(self.yc[np.asarray(self.wing_jy)], dtype=float)
        xrel = np.asarray(self.x_ctrl, dtype=float)      # 0..1 local chord
        nc, ng = len(xrel), len(ys)
        if nc < 2 or ng < 1:
            return None
        # ---- per-span station geometry (planform) ---------------------
        # xle_l / c_l / z_l per station; the stock case is (0, chord, 0),
        # which reproduces the single-section lattice exactly.
        xle_l = np.zeros(ng)
        c_l = np.full(ng, float(self.chord))
        z_l = np.zeros(ng)
        pl = self.planform
        if pl is not None and 'c_j' in pl:
            _jy = np.asarray(self.wing_jy)
            xle_l = np.asarray(pl['xle_j'], float)[_jy]
            c_l = np.maximum(np.asarray(pl['c_j'], float)[_jy], 1e-9)
            z_l = np.asarray(pl['z_j'], float)[_jy]
        xs = xrel * float(self.chord)
        try:
            nwk = max(1, int(os.environ.get('WING3D_KUTTA_KERNEL_WKP', '4')))
            wlen = float(os.environ.get('WING3D_KUTTA_KERNEL_WLEN', '3'))
            nsplit = max(1, int(os.environ.get(
                'WING3D_KUTTA_KERNEL_SPLIT', '1')))
        except ValueError:
            nwk, wlen, nsplit = 4, 3.0, 1
        # cells: chordwise between STATIONS (cell i carries the loading
        # Gamma_i = mu_i - mu_{i-1}, so its chord extent is exactly
        # [x_{i-1}, x_i], x_{-1} = LE), spanwise centred on a station.
        exr = np.empty(nc + 1)      # RELATIVE cell edges (0 .. 1)
        exr[0] = 0.0
        exr[1:] = xrel
        ey = np.empty(ng + 1)
        ey[0] = -0.5 * self.span
        ey[1:ng] = 0.5 * (ys[:-1] + ys[1:])
        ey[ng] = 0.5 * self.span
        exa = xle_l[:, None] + c_l[:, None] * exr[None, :]   # absolute x
        # horseshoe per cell: bound segment at the 1/4 chord, two legs to a
        # far downstream station (the trailed wake is IN the kernel, so its
        # downwash is inside the tangency equation instead of being smeared
        # into the field).
        fars = [self.chord]
        dxw = 0.5 * wlen * self.span / nwk
        for _ in range(nwk):
            fars.append(fars[-1] + dxw)
            dxw *= 1.6
        seg = []          # (nseg, 2, 3) filament endpoints
        col = []          # (seg,) station cell index the filament belongs to
        wgt = []          # (seg,) share of that cell's loading it carries
        xcp = np.empty(nc)          # root-frame chordwise control point
        xcpr = np.empty((nc, ng))   # 3/4 chord in LOCAL chord units
        ycp = 0.5 * (ey[:-1] + ey[1:])   # spanwise control point (mid-cell)
        zcp = z_l.copy()
        for i in range(nc):
            xcp[i] = self.chord * (0.25 * exr[i] + 0.75 * exr[i + 1])
            for l in range(ng):
                y0, y1 = ey[l], ey[l + 1]
                z0 = z_l[l]
                z1 = z_l[min(l + 1, ng - 1)]
                l1 = min(l + 1, ng - 1)
                # chordwise cell [exr[i], exr[i+1]] at station l and at l+1
                xa0, xb0 = exa[l, i], exa[l, i + 1]
                xa1, xb1 = exa[l1, i], exa[l1, i + 1]
                xcpr[i, l] = 0.25 * exr[i] + 0.75 * exr[i + 1]
                c = i * ng + l
                # The station cell carries the loading Gamma_i = mu_i -
                # mu_{i-1} spread over [x_{i-1}, x_i]; distributing it over
                # WING3D_KUTTA_KERNEL_SPLIT sub-horseshoes is a refinement
                # of the same density, not a change of unknown.
                for k in range(nsplit):
                    # sub-cell of the chordwise cell, same at both stations
                    xqs0 = (xa0 + (xb0 - xa0) * k / nsplit
                            + 0.25 * (xb0 - xa0) / nsplit)
                    xqs1 = (xa1 + (xb1 - xa1) * k / nsplit
                            + 0.25 * (xb1 - xa1) / nsplit)
                    f = 1.0 / nsplit
                    # bound vortex, +y (so +Gamma gives downwash aft of it)
                    seg.append([[xqs0, y0, z0], [xqs1, y1, z1]])
                    col.append(c)
                    wgt.append(f)
                    # trailing legs: a CLOSED horseshoe, so the y0 leg runs
                    # upstream (-x) and the y1 leg downstream (+x). Adjacent
                    # span cells then leave (Gamma_{l-1}-Gamma_l) on their
                    # shared leg -- the classical VLM wake -- instead of
                    # summing (a same-direction bug measured as -0.76*alpha
                    # instead of -alpha on the exact flat plate).
                    for yy, zz, xq, flip in ((y0, z0, xqs0, True),
                                             (y1, z1, xqs1, False)):
                        for kk in range(len(fars) - 1):
                            a0 = xq if kk == 0 else fars[kk]
                            seg.append([[fars[kk + 1], yy, zz],
                                        [a0, yy, zz]] if flip else
                                       [[a0, yy, zz],
                                        [fars[kk + 1], yy, zz]])
                            col.append(c)
                            wgt.append(f)
        info = dict(nc=nc, ng=ng, seg=np.asarray(seg), col=np.asarray(col),
                    wgt=np.asarray(wgt), exr=exr, ey=ey, xs=xs, ys=ys,
                    xcp=xcp, xcpr=xcpr, ycp=ycp, zcp=zcp, xf=fars[-1],
                    nwk=nwk, wlen=wlen, nsplit=nsplit)
        return info

    def dlm_kernel(self, rebuild=False):
        """Analytic doublet-lattice downwash matrix D (mode 6).

        D[(i,j),(k,l)] = normal (z) velocity at the control point of cell
        (i,j) induced by a UNIT loading of cell (k,l), i.e. by the horseshoe
        a constant doublet cell is equivalent to. Shape (nc*ng, nc*ng),
        i-major to match mu_b.ravel().

        The cell LOADING is the primary quantity in this form:

            Gamma_cell(i,l) = mu(i,l) - mu(i-1,l),    mu(-1) = 0

        so mu(TE) is the cumulative sum of the cell loadings and there is no
        potential jump at the leading edge BY CONSTRUCTION -- exactly the
        condition the exact flat plate satisfies (mu(0) = 0). That is what
        the constant-doublet form (mu primary, one constant panel reaching
        x=0) could not do: its first panel dumps the whole mu(0.043c) jump
        on the LE edge, and the LE tangency row then reads +0.41*alpha where
        the exact solution needs -alpha (measured, tmp_a12_dbg2).

        SIGN: the bound filament runs +y and the legs trail downstream, so a
        positive loading gives DOWNWASH aft of it (physical lift). Verified
        against the exact 2D flat plate in tests3d/tsd_a12_kernel.py: with
        mu = pi*alpha*F(x) the kernel returns w = -alpha along the chord.
        """
        if getattr(self, '_dlmK', None) is not None and not rebuild:
            return self._dlmK
        info = self._dlm_panels()
        self._dlmInfo = info
        if info is None:
            self._dlmK = np.zeros((0, 0))
            return self._dlmK
        from .ckernels import loop_vel_block
        nc, ng = info['nc'], info['ng']
        seg = info['seg']
        ed = seg[:, None, :, :]        # (nseg, 1, 2, 3)
        emask = np.ones((len(seg), 1), dtype=bool)
        P = np.zeros((nc * ng, 3))
        _pl = self.planform if (self.planform is not None
                                and 'c_j' in self.planform) else None
        for i in range(nc):
            for j in range(ng):
                if _pl is None:
                    P[i * ng + j] = (info['xcp'][i], info['ycp'][j], 0.0)
                else:
                    jj = int(np.asarray(self.wing_jy)[j])
                    P[i * ng + j] = (
                        float(_pl['xle_j'][jj])
                        + info['xcpr'][i, j] * float(_pl['c_j'][jj]),
                        info['ycp'][j], info['zcp'][j])
        if os.environ.get('WING3D_KUTTA_KERNEL_PG', '0') == '1':
            # Prandtl-Glauert stretch (x -> x/beta) as an ABLATION: the
            # kernel is geometric by default so every Mach effect in the
            # answer comes from the field, not from a linearised similarity
            # that is singular at M=1.
            beta = max(np.sqrt(max(1.0 - self.Minf ** 2, 1e-6)), 0.3)
            P[:, 0] /= beta
            ed = ed.copy()
            ed[..., 0] /= beta
        w = loop_vel_block(P, ed, emask)[:, :, 2]     # (ncp, nseg)
        S = np.zeros((len(seg), nc * ng))
        S[np.arange(len(seg)), info['col']] = info['wgt']
        DG = w @ S                                    # (ncp, ncell)
        # loading -> doublet: Gamma(i,l) = mu(i,l) - mu(i-1,l)
        C = np.eye(nc * ng) - np.eye(nc * ng, k=-ng)
        self._dlmK = DG @ C
        return self._dlmK

    def dlm_tangency_residual(self):
        """Mode 6 closure: the tangency equation, LINEAR in mu.

            D . mu = dyc/dx(x_cp_i) - alpha   at the 3/4 chord of each cell

        Same sign convention as the mode-5 row (`wav - (slope_c - alpha)`):
        the sheet's analytic downwash must equal the geometric slope minus
        the angle of attack (alpha>0 -> downwash -> D.mu<0). Row nc-1 is the
        discrete Kutta condition in mu form, mu(TE) = mu(last-1), unchanged
        from A11. Returns (nc, ng), scaled by kutta_w in the coupled
        residual.
        """
        nc = self.n_ctrl
        jy = np.asarray(self.wing_jy)
        ng = len(jy)
        if nc == 0 or ng == 0 or self.mu_b.size == 0:
            return np.zeros((nc, ng))
        D = self.dlm_kernel()
        info = self._dlmInfo
        slope_c = self._camber_at(np.asarray(info['xcpr'], float))
        g = self._geom_target(slope_c, np.asarray(info['xcp'], float))
        w = (D @ self.mu_b.ravel()).reshape(nc, ng)
        if self.kutta_te_kutta:
            # Kutta row in mu units: mu(TE) - mu(last-1) = 0
            w[-1, :] = self.mu_b[-1, :] - self.mu_b[-2, :]
            g[-1, :] = 0.0
        w, g = self._zero_rows(w, g)
        return w - g

    # ------------------------------------------------------------------
    # mode 7: HYBRID row -- w_field(phi) - (E - D).mu
    # ------------------------------------------------------------------
    def _stencil_heights(self):
        """The two z that `w_plane` actually samples (its own geometry).

        wu = (phi[k0+1]-phi[k0])/(z[k0+1]-z[k0]) is a first-order
        difference about the cell centre (z[k0]+z[k0+1])/2, and wl about
        (z[k0-1]+z[k0-2])/2; a one-sided difference of a smooth potential
        IS the normal velocity at the cell centre. So those two heights are
        where the estimator reads w, and they are what E must be sampled
        at for the bias to cancel. On g80 they are +0.0074c (above the
        sheet) and -0.0246c (below), i.e. NOT symmetric: that asymmetry is
        the estimator's, and it has to be reproduced, not idealised.
        """
        k0 = self.k0
        kp = min(k0 + 1, self.nz - 1)
        km = max(k0 - 2, 0)
        return (0.5 * (self.z[kp] + self.z[k0]),
                0.5 * (self.z[k0 - 1] + self.z[km]))

    def dlm_field_equivalent(self, rebuild=False):
        """E_node = `w_plane`'s estimator applied to the analytic sheet.

        Rows are the FIELD's own wing nodes (x_c[iu], y_c[jy]) -- exactly
        the points `w_plane()` returns -- and columns are the same mu
        layout as D. Shape (n_node*ng, nc*ng). Each row is the mean of the
        horseshoe downwash sampled at `_stencil_heights()`, i.e. the same
        two-z average `wav = 0.5*(wu+wl)` takes, built from the same
        Biot-Savart source as D (`ckernels.loop_vel_block`). It is a pure
        function of the grid and the panel lattice: no solution enters.
        """
        if getattr(self, '_dlmE', None) is not None and not rebuild:
            return self._dlmE
        if getattr(self, '_dlmInfo', None) is None:
            self.dlm_kernel()
        info = self._dlmInfo
        if info is None:
            self._dlmE = np.zeros((0, 0))
            return self._dlmE
        from .ckernels import loop_vel_block
        nc, ng = info['nc'], info['ng']
        seg = info['seg']
        ed = seg[:, None, :, :]
        emask = np.ones((len(seg), 1), dtype=bool)
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        xN = np.asarray(self.xc[iu], dtype=float)
        yN = np.asarray(self.yc[jy], dtype=float)
        z_up, z_lo = self._stencil_heights()
        nn = len(xN)
        P = np.zeros((nn * ng, 3))
        for i in range(nn):
            P[i * ng:(i + 1) * ng, 0] = xN[i]
            P[i * ng:(i + 1) * ng, 1] = yN
        P[:, 2] = z_up
        w = loop_vel_block(P, ed, emask)[:, :, 2]
        P[:, 2] = z_lo
        w = w + loop_vel_block(P, ed, emask)[:, :, 2]
        S = np.zeros((len(seg), nc * ng))
        S[np.arange(len(seg)), info['col']] = info['wgt']
        DG = (0.5 * w) @ S
        C = np.eye(nc * ng) - np.eye(nc * ng, k=-ng)
        self._dlmE = DG @ C
        return self._dlmE

    def set_thickness_reference(self, w_node):
        """Attach the circulation-OFF reference field (A14, B' = B + w_ref).

        `w_node` is `w_plane`'s mean-plane normal velocity (wav) on the
        FIELD's own wing nodes, shape (n_node, ng), measured on a solve of
        the SAME discrete problem with the bound/wake doublet sheets
        disabled. In the continuum a thickness-only (source) transpiration
        has wav = 0 at the mean plane, so this vector is a pure
        discretisation/reference offset -- exactly the part of `w_field`
        that the doublet sheet does not explain and that the Mach signal
        was buried under (A13 §6). It is mu-INDEPENDENT, so folding it into
        the bias keeps the row linear in mu.
        """
        self._thickRefNodes = None if w_node is None \
            else np.asarray(w_node, dtype=float)

    def dlm_thickness_reference(self):
        """w_ref = L . (reference wav on the field nodes), or None."""
        if getattr(self, '_thickRefNodes', None) is None:
            return None
        if getattr(self, '_dlmInfo', None) is None:
            self.dlm_kernel()
        info = self._dlmInfo
        if info is None:
            return None
        nc, ng = info['nc'], info['ng']
        w = np.asarray(self._thickRefNodes, dtype=float)
        if w.shape != (len(self.wing_ix), ng):
            raise ValueError('thickness reference shape %s != %s' % (
                w.shape, (len(self.wing_ix), ng)))
        return (self.dlm_obs_matrix() @ w.ravel()).reshape(nc, ng)

    def dlm_obs_matrix(self, rebuild=False):
        """L: (nc*ng, n_node*ng) -- field node lattice -> lattice control pts.

        `w_plane` lives on the FIELD's wing nodes; the lattice collocates
        at the 3/4 chord of each station cell (D is validated there and
        must not move). L is the linear map between the two, and it is
        applied IDENTICALLY to the field estimator and to E, so it cancels
        exactly in B = L.E_node - D. Only the DIFFERENCE has to be smooth
        for the cancellation to hold, not either side on its own.
        Linear extrapolation is used past the ends (the LE cell control
        point 0.75*x_0 lies inboard of the first node and carries the
        singular loading, so clamping would be worse).
        """
        if getattr(self, '_dlmL', None) is not None and not rebuild:
            return self._dlmL
        if getattr(self, '_dlmInfo', None) is None:
            self.dlm_kernel()
        info = self._dlmInfo
        if info is None:
            self._dlmL = np.zeros((0, 0))
            return self._dlmL
        nc, ng = info['nc'], info['ng']
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        xN = np.asarray(self.xc[iu], dtype=float)
        yN = np.asarray(self.yc[jy], dtype=float)
        xq = np.asarray(info['xcp'], dtype=float)
        yq = np.asarray(info['ycp'], dtype=float)
        nn = len(xN)

        def lw(q, nodes):
            i1 = int(np.clip(np.searchsorted(nodes, q), 1, len(nodes) - 1))
            i0 = i1 - 1
            t = (q - nodes[i0]) / max(nodes[i1] - nodes[i0], 1e-30)
            return ((i0, 1.0 - t), (i1, t))

        L = np.zeros((nc * ng, nn * ng))
        for i in range(nc):
            (a0, wa0), (a1, wa1) = lw(xq[i], xN)
            for j in range(ng):
                (b0, wb0), (b1, wb1) = lw(yq[j], yN)
                r = i * ng + j
                for aa, wax in ((a0, wa0), (a1, wa1)):
                    for bb, wby in ((b0, wb0), (b1, wb1)):
                        L[r, aa * ng + bb] += wax * wby
        self._dlmL = L
        return L

    def dlm_hybrid_bias(self, rebuild=False):
        """B = L.E_node - D: the mode-5 estimator's bias as a FIXED matrix.

        B is grid- and kernel-determined (no solution enters), so the
        hybrid row

            w_field(phi) - B.mu - (dyc/dx - alpha) = 0

        is LINEAR in mu and depends on phi only through w_field -- the
        structure mode 6 lacked. In the limit where the field really is
        the analytic doublet sheet, w_field = L.E_node.mu and the row
        collapses to D.mu = dyc/dx - alpha exactly (mode 6, validated), so
        the bias cancels identically and any deviation is a physical
        difference the field carries (compressibility, thickness,
        farfield), not estimator bias.
        """
        if getattr(self, '_hybB', None) is not None and not rebuild:
            return self._hybB
        self._hybB = self.dlm_obs_matrix() @ self.dlm_field_equivalent() \
            - self.dlm_kernel()
        return self._hybB

    def hybrid_tangency_residual(self):
        """Mode 7 closure: bias-free chordwise shape AND a live phi loop.

            row[i,j] = wav_field(phi; x_cp_i, y_cp_j)
                       - (B.mu)[i,j] - (dyc/dx(x_cp_i) - alpha)

        The field-induced normal velocity is the PHYSICS (so Mach, shocks,
        entropy and the BL displacement all reach the circulation through
        w_field), while B removes the estimator's systematic chordwise
        bias against the validated kernel. Row nc-1 stays the discrete
        Kutta condition mu(TE) = mu(last-1), as in modes 5/6. Returns
        (nc, ng), scaled by kutta_w in the coupled residual.
        """
        nc = self.n_ctrl
        jy = np.asarray(self.wing_jy)
        ng = len(jy)
        if nc == 0 or ng == 0 or self.mu_b.size == 0:
            return np.zeros((nc, ng))
        B = self.dlm_hybrid_bias()
        info = self._dlmInfo
        slope_c = self._camber_at(np.asarray(info['xcpr'], float))
        g = self._geom_target(slope_c, np.asarray(info['xcp'], float))
        _wu, _wl, wav = self.w_plane()
        iu = np.asarray(self.wing_ix)
        w_node = np.asarray(wav[np.ix_(iu, jy)], dtype=float)
        w_field = (self.dlm_obs_matrix() @ w_node.ravel()).reshape(nc, ng)
        w = w_field - (B @ self.mu_b.ravel()).reshape(nc, ng)
        wref = self.dlm_thickness_reference()
        if wref is not None:
            # A14: B' = B + w_ref -- remove the circulation-off reference
            # field's mean-plane downwash (thickness response + box/farfield
            # + discrete-Green's-function offset), so the only phi-dependence
            # left in the row is the CIRCULATORY part of the field.
            w = w - wref
        if self.kutta_te_kutta:
            w[-1, :] = self.mu_b[-1, :] - self.mu_b[-2, :]
            g[-1, :] = 0.0
        w, g = self._zero_rows(w, g)
        return w - g

    def hybrid_coupling_report(self, eps=1e-7):
        """A13 diagnostic: is the mu block really phi-coupled?

        Measures, on the ASSEMBLED rows (not on the operators):
          dphi     = max |d(row)/d(phi_node)| by central difference on the
                     four stencil nodes of one control point, and
          dmu      = max |d(row)/d(mu)| by central difference on one
                     diagonal mu DOF,
        plus the analytic 1-norms of the corresponding matrix blocks. The
        mode-6 comparison calls dlm_tangency_residual(), whose phi block is
        structurally zero (that is the defect mode 7 removes).
        """
        nc, ng = self.n_ctrl, len(self.wing_jy)
        out = {'nc': int(nc), 'ng': int(ng)}
        if nc == 0 or ng == 0 or self.mu_b.size == 0:
            return out
        B = self.dlm_hybrid_bias()
        out['B_max'] = float(np.abs(B).max())
        out['B_l1_row_max'] = float(np.abs(B).sum(axis=1).max())
        out['D_max'] = float(np.abs(self.dlm_kernel()).max())
        ic, jc = nc // 2, ng // 2
        k0 = self.k0
        nodes = [(min(k0 + 1, self.nz - 1), 1.0), (k0, -1.0),
                 (k0 - 1, 1.0), (max(k0 - 2, 0), -1.0)]
        iu = int(np.asarray(self.wing_ix)[self.i_ctrl[ic]])
        jy = int(np.asarray(self.wing_jy)[jc])
        base = None
        best = 0.0
        for kk, _s in nodes:
            old = float(self.phi[iu, jy, kk])
            self.phi[iu, jy, kk] = old + eps
            r1 = float(self.hybrid_tangency_residual()[ic, jc])
            self.phi[iu, jy, kk] = old - eps
            r2 = float(self.hybrid_tangency_residual()[ic, jc])
            self.phi[iu, jy, kk] = old
            if base is None:
                base = 0.5 * (r1 + r2)
            best = max(best, abs(r1 - r2) / (2.0 * eps))
        out['dphi_row'] = float(best)
        r0 = self.hybrid_tangency_residual()[ic, jc]
        mu0 = float(self.mu_b[ic, jc])
        self.mu_b[ic, jc] = mu0 + eps
        r1 = float(self.hybrid_tangency_residual()[ic, jc])
        self.mu_b[ic, jc] = mu0 - eps
        r2 = float(self.hybrid_tangency_residual()[ic, jc])
        self.mu_b[ic, jc] = mu0
        out['dmu_row'] = float(abs(r1 - r2) / (2.0 * eps))
        out['row'] = float(r0)
        # mode-6 control: this must be identically zero
        zp = np.zeros(4)
        for n, (kk, _s) in enumerate(nodes):
            old = float(self.phi[iu, jy, kk])
            self.phi[iu, jy, kk] = old + eps
            zp[n] = float(self.dlm_tangency_residual()[ic, jc])
            self.phi[iu, jy, kk] = old
        out['dphi_row_mode6'] = float(abs(zp - zp[0]).max() / eps)
        out['ratio_phi_mu'] = out['dphi_row'] / max(out['dmu_row'], 1e-30)
        return out

    def tangency_panel_residual(self):
        """Mode 5 closure: one tangency row per (panel, span station).

            row[i, j] = wav(x_i, y_j) - (dyc/dx(x_i) - alpha)

        i runs over the nc chordwise CONTROL stations (a subset of the
        field's own wing nodes), j over the span cells. This is the
        chordwise tangency problem itself, not a reduction of it: every
        chordwise panel sees its own local surface slope, so the LE
        region and the aft loading are both represented and the shape of
        the chordwise loading is the solution instead of a prescribed
        F(x). The TE row asks for zero aft loading at the last station --
        that IS the Kutta condition in doublet form (dmu/dx -> 0 at the
        TE).

        Returns shape (nc, ng), scaled by kutta_w in the coupled
        residual; flattened (i-major) to match the mu_b unknown layout.
        """
        k0 = self.k0
        iu = np.asarray(self.wing_ix)
        jy = np.asarray(self.wing_jy)
        if len(iu) == 0 or len(jy) == 0 or self.n_ctrl == 0:
            return np.zeros((0, len(jy)))
        xc = self.xc[iu] / self.chord
        _wu, _wl, wav = self.w_plane()
        rows = np.zeros((self.n_ctrl, len(jy)))
        for p in range(self.n_ctrl):
            ip = int(iu[self.i_ctrl[p]])
            x_at = float(xc[self.i_ctrl[p]])
            if self.profile is not None:
                slope_c = float(self.profile['camber'](x_at))
            else:
                slope_c = 0.0
            rows[p, :] = wav[ip, jy] - self._geom_target(
                np.array([slope_c]), np.array([float(self.xc[ip])]))[0]
        return rows

    def mu_circulation(self):
        """Section circulation Gamma(y) implied by the mu matrix.

        FOR A CONTINUOUS doublet sheet the loading is gamma = dmu/dx and
        the section circulation is gamma(1) - gamma(0) = mu(TE) - mu(LE);
        with a discrete mu the same statement is the sum of the
        chordwise increments. Reported so the mode-5 circulation can be
        compared with the lifting line without a weight (CL_circ =
        2*sum(Gamma*dy)/S, the same definition mode 4 uses).
        """
        if self.mu_b.size == 0:
            return np.zeros(len(self.wing_jy))
        if self.kutta_kernel:
            # mode 6: the cell loading is Gamma_i = mu_i - mu_{i-1} with
            # mu(x=0) = 0, so the section circulation is mu(TE) itself.
            return self.mu_b[-1, :].copy()
        return self.mu_b[-1, :] - self.mu_b[0, :]

    def mu_loading(self):
        """Chordwise loading gamma(x) = dmu/dx at each span station.

        Returns (nc, ng) per-unit-length loading from the discrete mu
        (forward differences between control stations). This is the
        quantity a physical Kutta-satisfying section must show: peak
        near the LE, falling to ~0 at the TE.
        """
        if self.mu_b.size == 0:
            return np.zeros((0, len(self.wing_jy)))
        xc = np.asarray(self.x_ctrl)
        out = np.zeros_like(self.mu_b)
        if len(xc) > 1:
            out[:-1, :] = (np.diff(self.mu_b, axis=0)
                           / np.diff(xc)[:, None])
            out[-1, :] = out[-2, :]
        return out

    def recirc_update(self, W=0.25, circff=None):
        """TSFOIL-style Kutta update (RECIRC).

        TE jump CIRCTE(y) from potentials EXTRAPOLATED to the TE point
        from both sides (CJ-style linear extrapolation -- avoids the
        cut singularity that poisons raw discrete jumps). Relaxed:
        Gamma <- (1-W)*Gamma + W*CIRCTE. Wake jump interpolated
        CIRCTE -> CIRCFF downstream (box consistency). Returns max |dGamma|.
        """
        k0 = self.k0
        # linear extrapolation weights to z=0 from planes k0+1/+2, k0-1/-2
        zp1, zp2 = self.z[min(k0 + 1, self.nz - 1)], self.z[min(k0 + 2, self.nz - 1)]
        zm1, zm2 = self.z[max(k0 - 1, 0)], self.z[max(k0 - 2, 0)]
        wu = abs(zp2) / max(abs(zp2 - zp1), 1e-12)
        wl = abs(zm2) / max(abs(zm2 - zm1), 1e-12)
        ite = self.wing_ix[-1]
        ku1 = min(k0 + 1, self.nz - 1)
        ku2 = min(k0 + 2, self.nz - 1)
        kl1 = max(k0 - 1, 0)
        kl2 = max(k0 - 2, 0)
        pup = wu * self.phi[ite, :, ku1] - (wu - 1.0) * self.phi[ite, :, ku2]
        plo = wl * self.phi[ite, :, kl1] - (wl - 1.0) * self.phi[ite, :, kl2]
        circte = pup - plo
        # map TE column (full ny) onto wake span stations
        jy = np.array(self.wing_jy)
        newG = np.interp(self.yc, self.yc[jy], circte[jy],
                         left=0.0, right=0.0)
        oldG = self.wake_G.copy()
        self.wake_G = (1 - W) * oldG + W * newG
        if circff is not None:
            self._circff = np.asarray(circff, dtype=float)
        return float(np.abs(self.wake_G - oldG).max())

    def line_relax(self, omega=0.8, direction=+1):
        """One x-line-relaxation sweep (all y,z lines, Thomas solve).

        Approximate-Newton step with the |Jacobian| tridiagonal in x
        (absolute values -> M-matrix, Thomas-stable even where the
        frozen A < 0; y/z/transpiration/cut explicit in RHS).
        Matches the SOR sign convention (x += +w*R/D).
        `direction` alternates sweep order (helps supersonic advection).
        Returns max |update|.
        """
        phi, nx, ny, nz = self.phi, self.nx, self.ny, self.nz
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
        phix_f = (phi[1:, :, :] - phi[:-1, :, :]) / dx[:, None, None]
        if self.linear:
            A_f = np.full_like(phix_f, max(1 - self.Minf ** 2, 1e-6))
        elif self._frozenA is not None:
            A_f = self._frozenA
        else:
            A_f = self._face_coeff()
        maxupd = 0.0
        R = self.residual()
        # Jacobi-style: buffer all line updates from ONE residual, apply
        # at sweep end (Gauss-Seidel with per-line stale R diverges).
        dall = np.zeros_like(phi)
        jr = range(1, ny - 1) if direction > 0 else range(ny - 2, 0, -1)
        kr = range(1, nz - 1) if direction > 0 else range(nz - 2, 0, -1)
        for j in jr:
            dyy = (1.0 / (dy[j - 1] * wy[j]) + 1.0 / (dy[j] * wy[j])) \
                if 0 < j < ny - 1 else 0.0
            for k in kr:
                dzz = (1.0 / (dz[k - 1] * wz[k]) + 1.0 / (dz[k] * wz[k])) \
                    if 0 < k < nz - 1 else 0.0
                n = nx - 2
                a = np.zeros(n)
                b = np.zeros(n)
                c = np.zeros(n)
                rhs = np.zeros(n)
                for ii, i in enumerate(range(1, nx - 1)):
                    wxi = wx[i]
                    A0 = abs(A_f[i - 1, j, k]) if i - 1 < A_f.shape[0] \
                        else 0.0
                    A1 = abs(A_f[i, j, k]) if i < A_f.shape[0] else 0.0
                    d0 = max(dx[i - 1], 1e-300)
                    d1 = max(dx[i] if i < len(dx) else dx[-1], 1e-300)
                    if ii > 0:
                        a[ii] = -A0 / (d0 * wxi)
                    b[ii] = (A0 / (d0 * wxi) + A1 / (d1 * wxi)
                             + dyy + dzz)
                    if ii < n - 1:
                        c[ii] = -A1 / (d1 * wxi)
                    rhs[ii] = R[i, j, k]
                for ii in range(1, n):
                    w = a[ii] / max(b[ii - 1], 1e-300)
                    b[ii] -= w * c[ii - 1]
                    rhs[ii] -= w * rhs[ii - 1]
                dphi = np.zeros(n)
                dphi[-1] = rhs[-1] / max(b[-1], 1e-300)
                for ii in range(n - 2, -1, -1):
                    dphi[ii] = (rhs[ii] - c[ii] * dphi[ii + 1]) / max(
                        b[ii], 1e-300)
                # match SOR sign convention (x += +w*R/D, A_ii<0)
                upd = omega * dphi
                for ii, i in enumerate(range(1, nx - 1)):
                    dall[i, j, k] = upd[ii]
                    if abs(upd[ii]) > maxupd:
                        maxupd = abs(upd[ii])
        self._apply_bc(dall)
        phi += dall
        return maxupd


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
