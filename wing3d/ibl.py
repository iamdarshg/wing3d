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
    # Re_theta at transition (Michel); guard Re_x=0 (stagnation) -> huge
    # threshold (no immediate transition)
    Rx = max(Re_x, 1.0)
    return 1.174 * (1.0 + 22400.0 / Rx ** 0.625) * Rx ** 0.46


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


def solve_streamline(s, Ue, nu, Re_ref=1.0, x_tr_forced=None):
    """Integrate IBL along arc-length stations s with edge velocity Ue.

    s, Ue: 1D arrays (same length), Ue > 0. nu: kinematic viscosity in the
    same units (Vinf = 1). Returns dict of arrays + transition/separation idx.
    """
    s = np.asarray(s, dtype=float)
    Ue = np.maximum(np.asarray(Ue, dtype=float), 1e-9)
    n = len(s)
    ds = np.diff(s)
    ds = np.append(ds, ds[-1] if len(ds) else 1.0)
    dUe = np.gradient(Ue, s)
    Re_x = np.cumsum(np.abs(np.diff(s, prepend=0)) * Ue) / max(nu, 1e-300)

    theta = np.zeros(n)
    H = np.zeros(n)
    cf = np.zeros(n)
    dstar = np.zeros(n)
    lam_arr = np.zeros(n)
    state = np.zeros(n, dtype=int)  # 0 laminar, 1 turbulent
    i_tr = n  # transition index (default: never)
    i_sep = -1

    # laminar start (stagnation-like): theta^2 = 0.075*nu*ds/Ue
    t2 = 0.075 * nu * max(ds[0], 1e-12) / Ue[0]
    bubble = False
    for i in range(n):
        if i > 0:
            lam = t2 / max(nu, 1e-300) * dUe[max(i - 1, 0)]
            lam_arr[i] = lam
            if lam < -0.09:
                # laminar separation before transition: separation-induced
                # transition (short bubble model -> turbulent from here)
                bubble = True
                i_tr = i
                break
            # clamp to Thwaites validity range (explicit-Euler overshoot
            # near stagnation otherwise collapses theta -> Cf ~ 1000)
            lam_c = min(lam, 0.12)
            F = thwaites_rhs(lam_c)
            t2 = max(t2 + F * nu / max(Ue[i], 1e-9) * ds[i], 1e-18)
        th = np.sqrt(max(t2, 1e-18))
        theta[i] = th
        Re_t = Ue[i] * th / max(nu, 1e-300)
        ell = thwaites_cf_ell(min(lam_arr[i], 0.12))
        cf[i] = min(2 * nu * ell / max(Ue[i] * th, 1e-300), 0.02)
        H[i] = 2.2 if lam_arr[i] > -0.02 else 2.6
        dstar[i] = H[i] * th
        # Michel transition check
        trforced = (x_tr_forced is not None and s[i] >= x_tr_forced)
        if Re_t > michel_transition(Re_x[i]) or trforced:
            i_tr = i
            break
    else:
        return {'s': s, 'Ue': Ue, 'theta': theta, 'H': H, 'cf': cf,
                'dstar': dstar, 'state': state, 'i_tr': n, 'i_sep': -1,
                'separated': False, 'bubble': False}

    # turbulent (Head) from transition with initial H
    H1 = 5.0
    th = theta[max(i_tr - 1, 0)] or 1e-6
    for i in range(i_tr, n):
        state[i] = 1
        if i > i_tr:
            # entrainment ODE for (Ue*theta*H1)
            F = head_F(H1)
            y = Ue[i - 1] * th * H1 + F * Ue[i - 1] * ds[i]
            H1 = y / max(Ue[i] * th, 1e-300)
            H1 = min(max(H1, 3.01), 12.0)
        H[i] = head_H_of_H1(H1)
        Re_t = Ue[i] * th / max(nu, 1e-300)
        cf[i] = ludwieg_tillmann_cf(Re_t, H[i])
        # momentum integral update for theta
        Cf = cf[i]
        dth = (Cf / 2 - th * (H[i] + 2) / max(Ue[i], 1e-9) * dUe[i]) \
            / max(Ue[i], 1e-9) * ds[i]
        th = max(th + dth, 1e-9)
        theta[i] = th
        dstar[i] = H[i] * th
        if H[i] > 2.4 or Cf <= 0:
            i_sep = i
            break
    return {'s': s, 'Ue': Ue, 'theta': theta, 'H': H, 'cf': cf,
            'dstar': dstar, 'state': state, 'i_tr': i_tr, 'i_sep': i_sep,
            'separated': (i_sep >= 0) or bubble, 'bubble': bubble}
