"""3D integral boundary layer along surface streamlines.

Laminar: Thwaites. Transition: Michel. Turbulent: Head entrainment +
Ludwieg-Tillmann skin friction. Small-crossflow approximation (streamline
coordinates), standard first-order 3D BL approach (cf. Drela 3D IBL work,
Cebeci-Smith 3D). Transition/turbulence act along each surface streamline
traced from the 3D panel solution -- no 2D strip solver involved.
"""
import numpy as np


def thwaites_rhs(lam):
    # Thwaites F(lambda); separation at lambda < -0.09
    return 0.45 - 6.0 * lam


def thwaites_cf_ell(lam):
    # l(lambda)*Re_theta approximation (Cebeci-Smith fit)
    return 0.22 + 1.57 * lam - 1.8 * lam ** 2


def michel_transition(Re_x):
    # Re_theta at transition (Michel 1951): 1.174*(1 + 22400/Re_x)*Re_x^0.46.
    # NOTE (fix): the 22400 term has NO ^0.625 exponent. With the exponent
    # the threshold sits ~7x too high (flat-plate transition at 1.3e8!);
    # without it, flat-plate transition lands at Re_x ~ 2e6 (textbook)
    # and airfoil transition moves with Re like TSFOIL2 (6-56%).
    # Guard Re_x=0 (stagnation) -> huge threshold (no immediate trip).
    Rx = max(Re_x, 1.0)
    return 1.174 * (1.0 + 22400.0 / Rx) * Rx ** 0.46


def head_F(H1):
    return 0.0306 * max(H1 - 3.0, 0.01) ** -0.573


def head_H_of_H1(H1):
    # H(H1) closure (Head/entrainment family fit)
    if H1 < 3.3:
        return 1.0 + 1.0 / max(H1 - 1.5, 0.2)
    return 1.0 + 1.5 / max(H1 - 2.0, 0.2)


def ludwieg_tillmann_cf(Re_theta, H):
    cf = 0.246 * 10.0 ** (-0.678 * H) * max(Re_theta, 1.0) ** -0.268
    return max(cf, 0.0)


def ludwieg_tillmann_cf(Re_theta, H):
    cf = 0.246 * 10.0 ** (-0.678 * H) * max(Re_theta, 10.0) ** -0.268
    return max(cf, 0.0)


def white_lambda_corr(lam):
    """White (2006) polynomial fits to Thwaites' table.

    l(lambda): dimensionless wall shear, cf = 2*nu*l/(Ue*theta).
    H(lambda): shape factor. Valid -0.09 <= lam <= 0.25.
    (Replaces Cebeci-Smith fits; continuous H vs old 2.2/2.6 step.)
    """
    l = 0.22 + 1.402 * lam + 0.018 * lam / (0.107 + lam)
    z = 0.25 - lam
    H = (2.0 + 4.14 * z - 83.5 * z ** 2 + 854.0 * z ** 3
         - 3337.0 * z ** 4 + 4576.0 * z ** 5)
    return l, min(max(H, 1.0), 10.0)


def head_H_to_H1(H):
    """Head (1958) exact H -> H1 (cf. TSFOIL2 IBL)."""
    H = max(H, 1.05)
    if H < 1.6:
        return 3.3 + 0.8234 * (H - 1.1) ** (-1.287)
    return 3.3 + 1.5501 * (H - 0.6778) ** (-3.064)


def head_H1_to_H(H1):
    """Head (1958) exact H1 -> H, two-branch."""
    H1 = max(H1, 3.01)
    dH1 = max(H1 - 3.3, 1e-6)
    Hb1 = 1.1 + (0.8234 / dH1) ** (1.0 / 1.287)
    if Hb1 <= 1.6:
        return max(Hb1, 1.05)
    return max(0.6778 + (1.5501 / dH1) ** (1.0 / 3.064), 1.05)


def solve_streamline(s, Ue, nu, Re_ref=1.0, x_tr_forced=None, me=None):
    """Integrate IBL along arc-length stations s with edge velocity Ue.

    s, Ue: 1D arrays (same length), Ue > 0. nu: kinematic viscosity in the
    same units (Vinf = 1). me: edge Mach (optional, for compressible
    momentum integral; defaults to 0 = incompressible). Returns dict of
    arrays + transition/separation idx.

    Laminar: Thwaites EXACT integral form (theta^2 = 0.45*nu/Ue^6 *
    int(Ue^5 ds)) with White (2006) l/H correlations -- no Euler
    marching, no overshoot collapse. Turbulent: Head entrainment via
    RK45 with exact H1/H closures (H0 = 1.4 canonical start).
    Separation bubbles (Horton): dead-air Cf=0 zone + theta doubling
    at reattachment.
    """
    s = np.asarray(s, dtype=float)
    Ue = np.maximum(np.asarray(Ue, dtype=float), 1e-9)
    n = len(s)
    if me is None:
        me = np.zeros(n)
    else:
        me = np.asarray(me, dtype=float)
    dUe = np.gradient(Ue, s)
    Re_x = np.cumsum(np.abs(np.diff(s, prepend=0)) * Ue) / max(nu, 1e-300)

    # --- laminar (Thwaites integral, exact) ---
    integ = np.zeros(n)
    ue5 = Ue ** 5
    for i in range(1, n):
        ds = max(s[i] - s[i - 1], 0.0)
        integ[i] = integ[i - 1] + 0.5 * (ue5[i] + ue5[i - 1]) * ds
    theta = np.sqrt(np.maximum(0.45 * nu * integ / Ue ** 6, 1e-30))
    lam_raw = theta ** 2 / max(nu, 1e-300) * dUe
    lam = np.clip(lam_raw, -0.09, 0.25)
    H = np.zeros(n)
    cf = np.zeros(n)
    for i in range(n):
        ell, H[i] = white_lambda_corr(lam[i])
        cf[i] = 2 * nu * ell / max(Ue[i] * theta[i], 1e-30)
    cf[0] = cf[1] if n > 1 else cf[0]  # remove stagnation singularity
    dstar = H * theta
    state = np.zeros(n, dtype=int)
    i_tr = n
    i_sep = -1
    bubble = False
    # laminar separation (unclipped lambda) -> Horton bubble
    sep_at = np.where(lam_raw < -0.09)[0]
    # Michel transition (skip LE singularity zone like TSFOIL2)
    Re_t = Ue * theta / max(nu, 1e-300)
    mic = np.array([michel_transition(rx) for rx in Re_x])
    trips = np.where(Re_t > mic)[0]
    trips = trips[trips >= 5]
    i_mich = int(trips[0]) if len(trips) else n
    if x_tr_forced is not None:
        cand = np.where(s >= x_tr_forced)[0]
        if len(cand):
            i_mich = min(i_mich, int(cand[0]))
    i_sep_cand = int(sep_at[0]) if len(sep_at) else n
    if i_sep_cand < i_mich:
        # separation-induced transition: Horton short bubble
        bubble = True
        i_tr = i_sep_cand
        i_sep = i_sep_cand
    else:
        i_tr = i_mich
    if i_tr >= n:
        return {'s': s, 'Ue': Ue, 'theta': theta, 'H': H, 'cf': cf,
                'dstar': dstar, 'state': state, 'i_tr': n, 'i_sep': -1,
                'separated': False, 'bubble': False}
    # Horton bubble zone (if separated): dead air then reattachment jump
    Lb = 0
    th_re = theta[max(i_tr - 1, 0)]
    H0 = 1.4
    if bubble:
        ue_sep = max(Ue[i_sep], 1e-9)
        # Horton short-bubble length scale (~5e4 * nu/Ue; capped)
        Lb = min(5e4 * nu / ue_sep, 0.3 * max(s[-1] - s[i_sep], 0.0))
        # dead-air zone: Cf ~ 0, frozen shape; reattachment doubles theta
        # (Horton: shear-layer entrainment across short bubbles)
        i_re = int(np.searchsorted(s, s[i_sep] + Lb))
        i_re = min(max(i_re, i_sep + 1), n - 1)
        cf[i_sep:i_re + 1] = 0.0
        th_re = 2.0 * max(theta[i_sep], 1e-10)
        H0 = 2.2  # elevated reattachment shape factor
        H[i_sep:i_re + 1] = H[i_sep]
        dstar[i_sep:i_re + 1] = H[i_sep] * theta[i_sep]
        i_tr = i_re
    # --- turbulent (Head entrainment, RK45) ---
    from scipy.integrate import solve_ivp

    def rhs(t, y):
        th, H1th = y
        th = max(th, 1e-10)
        H1 = max(H1th / th, 3.01)
        Hc = head_H1_to_H(H1)
        u_ = float(np.interp(t, s, Ue))
        du_ = float(np.interp(t, s, dUe))
        me_ = float(np.interp(t, s, me))
        Re_ = u_ * th / max(nu, 1e-300)
        cf_ = ludwieg_tillmann_cf(Re_, Hc)
        dth = 0.5 * cf_ - (Hc + 2.0 - me_ ** 2) * th / max(u_, 1e-9) * du_
        dH = 0.0306 * max(H1 - 3.0, 1e-9) ** (-0.6169)
        return [dth, dH]

    H1_0 = head_H_to_H1(H0)
    y0 = [max(th_re, 1e-10), H1_0 * max(th_re, 1e-10)]
    sol = solve_ivp(rhs, [s[i_tr], s[-1]], y0, method='RK45',
                    t_eval=s[i_tr:], rtol=1e-5, atol=1e-9)
    th_t = np.maximum(sol.y[0], 1e-10)
    H1_t = np.maximum(sol.y[1] / th_t, 3.01)
    H_t = np.array([head_H1_to_H(v) for v in H1_t])
    Re_t2 = Ue[i_tr:] * th_t / max(nu, 1e-300)
    cf_t = np.array([ludwieg_tillmann_cf(r_, h_)
                     for r_, h_ in zip(Re_t2, H_t)])
    theta[i_tr:] = th_t
    H[i_tr:] = H_t
    cf[i_tr:] = cf_t
    dstar[i_tr:] = H_t * th_t
    state[i_tr:] = 1
    # turbulent separation (H threshold / Cf floor)
    seps = np.where((H[i_tr:] > 2.4) | (cf[i_tr:] <= 0))[0]
    if len(seps):
        i_sep = int(i_tr + seps[0])
    return {'s': s, 'Ue': Ue, 'theta': theta, 'H': H, 'cf': cf,
            'dstar': dstar, 'state': state, 'i_tr': int(i_tr),
            'i_sep': int(i_sep),
            'separated': (i_sep >= 0) or bubble, 'bubble': bubble}
