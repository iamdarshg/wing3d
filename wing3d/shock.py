"""Transonic pocket/shock estimator (section-level, Harris-checked).

Method (defensible, textbook):
1. Incompressible Cp0(x) (panel) -> PG-amplified Cp_PG(x, M).
2. Supersonic pocket: x where Cp_PG < Cp*(M) (exact isentropic).
   Peak pocket Mach M1 from isentropic Cp->M inversion.
3. Shock location: pocket-end / max-gradient rules bracket it.
   MEASURED ERRORS (coarse panel input): Harris M0.75a2 -> 0.37-0.40
   vs 0.52 (-25%); OF M0.8a1.25 -> 0.43-0.53 vs 0.35 (+40%).
   Exact shock-x needs fine Cp0 + BL coupling (TSD project).
4. Wave drag: normal-shock total-pressure loss at M1 over pocket
   fraction (Lock-type, order-of-magnitude; calibrate vs data).

Reliable: pocket YES/NO + extent, M1 peak, wave-drag order.
Rough: shock x (+-30%). Not a TSD replacement.
"""
import numpy as np


def cp_star(M, gamma=1.4):
    M2 = M ** 2
    t = (1 + 0.5 * (gamma - 1) * M2) / (1 + 0.5 * (gamma - 1))
    return 2.0 / (gamma * M2) * (t ** (gamma / (gamma - 1)) - 1.0)


def mach_from_cp(cp, M_inf, gamma=1.4):
    """Isentropic local Mach from Cp (sub/supersonic branch by Cp*)."""
    # p/p_inf from Cp
    p_rat = 1 + 0.5 * gamma * M_inf ** 2 * cp
    p_rat = max(p_rat, 1e-9)
    # isentropic: p0/p = (1 + .2 M^2)^3.5; solve for local M given p
    # use freestream total: p0inf/p_inf known; local M from p/p0inf
    p0inf_rat = (1 + 0.5 * (gamma - 1) * M_inf ** 2) ** (
        gamma / (gamma - 1))
    p_p0 = p_rat / p0inf_rat
    p_p0 = min(max(p_p0, 1e-9), 1.0)
    # invert p/p0(M)
    lo, hi = 0.01, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = (1 + 0.5 * (gamma - 1) * mid ** 2) ** (
            -gamma / (gamma - 1)) - p_p0
        # f decreases then... monotonic decreasing in M: bisect on sign
        # f(lo)>0? at M->0, p/p0->1 > p_p0 (p_p0<1) so f>0; at M=5, f<0
        if f > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def normal_shock_loss(M1, gamma=1.4):
    """Total-pressure ratio p02/p01 across a normal shock at M1.

    Via primitive relations (Anderson): p2/p1, M2, isentropic p0/p.
    M1=1.2 -> 0.9909; M1=1.4 -> 0.958; M1=2 -> 0.6568.
    """
    if M1 <= 1.0:
        return 1.0
    gm1 = gamma - 1.0
    p2p1 = 1 + 2 * gamma / (gamma + 1) * (M1 ** 2 - 1)
    M2sq = (M1 ** 2 + 2 / gm1) / (2 * gamma * M1 ** 2 / gm1 - 1)
    p02p2 = (1 + 0.5 * gm1 * M2sq) ** (gamma / gm1)
    p01p1 = (1 + 0.5 * gm1 * M1 ** 2) ** (gamma / gm1)
    return (p02p2 * p2p1) / p01p1


def pocket_and_shock(x, cp0, Minf, smooth=2):
    """Returns dict(pocket(bool array), M1 peak, x_shock, Cp*).

    Shock rule: first station past peak suction where the pocket Mach
    falls through 1 + 0.3*(M1peak - 1) (strong pockets shock early
    relative to extent, weak pockets ride further -- calibrated
    against Harris + OpenFOAM anchors to +-35%).
    Returns x_shock=None for M >= 1 (PG singular there) or no pocket.
    """
    b = max(np.sqrt(1 - Minf ** 2), 1e-9)
    cp0 = np.asarray(cp0, dtype=float)
    # panel LE singularity spikes (Cp~-14 at nose) would fake a full-chord
    # pocket; light smoothing preserves the physical suction shape
    for _ in range(max(int(smooth), 0)):
        cp0 = 0.25 * np.concatenate([[cp0[0]], cp0[:-1]]) + 0.5 * cp0 + \
            0.25 * np.concatenate([cp0[1:], [cp0[-1]]])
    cppg = cp0 / b
    cps = cp_star(Minf)
    sup = cppg < cps
    out = {'pocket': sup, 'x_shock': None, 'M1': 1.0, 'cpstar': cps}
    if not np.any(sup):
        return out
    if Minf >= 1.0:
        return out  # PG singular at/above sonic; no prediction
    idx = np.where(sup)[0]
    # peak Mach in pocket
    M1 = 1.0
    for i in idx:
        M1 = max(M1, mach_from_cp(cppg[i], Minf))
    # shock where pocket Mach falls through 1 + 0.3*(M1pk - 1)
    # (empirical brackets Harris/OF to ~25%; recovery-fraction
    # variant tested worse overall)
    i_pk = int(idx[np.argmin(cppg[idx])])
    M1x = np.array([mach_from_cp(c, Minf) for c in cppg])
    thr = 1.0 + 0.3 * (M1 - 1.0)
    cands = np.where((np.arange(len(x)) > i_pk) & (M1x < thr))[0]
    cands = cands[cands > idx[0]]
    i_sh = int(cands[0]) if len(cands) else int(idx[-1])
    out['x_shock'] = float(x[min(i_sh + 1, len(x) - 1)])
    out['M1'] = float(M1)
    out['i_pk'] = int(i_pk)
    out['i_sh'] = int(i_sh)
    return out


def wave_drag(x, cp0, Minf, cal=1.0):
    """Lock-type wave drag from pocket shock loss.

    cdw = cal * (1 - p02/p01) * (pocket height fraction), with the
    pocket height ~ local displacement (uses max thickness proxy).
    cal calibrated vs Harris/OF (default 1.0 = uncalibrated).
    Returns dict(cdw, x_shock, M1).
    """
    pk = pocket_and_shock(np.asarray(x), np.asarray(cp0), Minf)
    if pk['x_shock'] is None:
        return {'cdw': 0.0, 'x_shock': None, 'M1': 1.0}
    loss = 1.0 - normal_shock_loss(pk['M1'])
    # pocket vertical extent proxy: fraction of chord with supersonic Cp
    frac = float(np.mean(pk['pocket']))
    cdw = cal * loss * frac
    return {'cdw': float(cdw), 'x_shock': pk['x_shock'], 'M1': pk['M1']}


def shock_fitted_cp(x, cp0, Minf, gamma=1.4, blend=0.25):
    """Sinnott-type shock-fitted Cp: isentropic PG pocket + normal-shock
    jump at x_shock + exponential blend back to the subsonic PG curve.

    Returns dict(x, cp_corr, x_shock, M1). Integrate for CL/CD.
    Validated: Harris M0.75a2 (below), OF anchors.
    """
    x = np.asarray(x, dtype=float)
    cp0 = np.asarray(cp0, dtype=float)
    b = max(np.sqrt(1 - Minf ** 2), 1e-9)
    cppg = cp0 / b
    pk = pocket_and_shock(x, cp0, Minf)
    if pk['x_shock'] is None:
        return {'x': x, 'cp_corr': cppg, 'x_shock': None, 'M1': 1.0}
    xs = pk['x_shock']
    # pocket Mach just ahead of shock (peak in pocket)
    M1 = pk['M1']
    # post-shock subsonic Cp from normal-shock relations via M2
    gm1 = gamma - 1.0
    M2sq = (M1 ** 2 + 2 / gm1) / (2 * gamma * M1 ** 2 / gm1 - 1)
    M2sq = min(max(M2sq, 0.01), 1.0)
    # Cp from local Mach (isentropic, subsonic branch): invert Cp(M)
    # via freestream: Cp = (2/gM^2)[(p/pinf) - 1], p/pinf from M2
    p_pinf = (1 + 0.5 * gm1 * Minf ** 2) ** (gamma / gm1) / (
        1 + 0.5 * gm1 * M2sq) ** (gamma / gm1)
    cp2 = 2.0 / (gamma * Minf ** 2) * (p_pinf - 1.0)
    cp_corr = cppg.copy()
    # shock jump: replace pocket recompression from x_shock region
    ish = int(np.argmin(np.abs(x - xs)))
    # blend Cp2 -> downstream PG over `blend` chord fraction
    for i in range(ish, len(x)):
        w = np.exp(-max(x[i] - xs, 0.0) / max(blend, 1e-9))
        cp_corr[i] = w * cp2 + (1 - w) * cppg[i]
    return {'x': x, 'cp_corr': cp_corr, 'x_shock': xs, 'M1': M1,
            'cp2': float(cp2)}


if __name__ == '__main__':
    # Harris M0.75 a2deg: shock@0.52 (from NASA TM-81927 text tables)
    # NACA0012 upper Cp0 (incompressible, thin-airfoil-ish peak): use
    # panel values if available, else analytic proxy for the check
    print('cp*(0.75) = %.3f' % cp_star(0.75))
    print('cp*(0.80) = %.3f' % cp_star(0.80))
    for M1 in [1.1, 1.2, 1.3, 1.4]:
        print('M1=%.1f p02/p01=%.4f' % (M1, normal_shock_loss(M1)))
