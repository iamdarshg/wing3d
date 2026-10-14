"""Transonic predictor for wing3d (empirical surrogate over OF anchors).

Covers the M0.65-1.1 gap where KT is invalid and TSD-nonlinear is
parked. Bilinear interpolation of the rhoCentralFoam Euler anchor
table (NACA0012 section, cases/openfoam/transonic/summary.json):
Mach nodes 0.8/0.9/0.95/1.0/1.1 at alpha=1.25deg, alpha nodes
0/1.25/4/8 at M0.8. Below M0.65 delegates to panel+KT (subsonic).

HONEST LIMITS (not physics): Euler (no BL decambering/friction --
add IBL CDf separately), NACA0012 section only (apply per-section
with sweep/3D corrections by caller), coarse-mesh anchors (no grid
bars; spurious ~0.04 drag floor at M0.8 included in CD -- subtract
for wave-drag estimates), alpha range 0-8deg, Mach 0.65-1.1.
Extrapolation outside returns NaN (no silent garbage).
"""
import numpy as np
import os
import json

_ANCHOR = {
    # (Mach, alpha): (CL, CD, shock_x, cpmin_up). Inflow-audited
    # 2026-10-08; N2-gas-consistent rows only. Finite-wing AR1.6
    # (LL-validated); Harris = independent 2D truth.
    (0.80, 0.00): (0.0003, 0.0474, None, None),
    (0.80, 1.25): (0.0686, 0.0485, 0.35, -1.27),
    (0.80, 4.00): (0.2187, 0.0464, 0.35, -0.96),
    (0.80, 8.00): (0.4040, 0.0415, 0.35, -1.27),
    (0.885, 1.25): (0.0696, 0.0512, 0.40, -0.74),
    (0.95, 1.25): (0.0830, 0.0650, 0.60, -0.76),
    (1.00, 1.25): (0.0790, 0.0770, 0.80, -0.72),
    (1.10, 1.25): (0.0533, 0.1030, 0.70, -0.60),
}
# Fine-mesh (86k-cell, AR1.6, a=1.25, rhoCentralFoam, N2 gas) campaign
# 2026-10-13: supersedes the coarse 12k-cell rows above. CL/CD here are
# the higher-fidelity anchors; CD floor ~0.038-0.040 vs old ~0.044
# spurious floor. Old rows kept for provenance/alpha-shape (a=0/4/8
# still only exist on the old M0.8 row).
_ANCHOR_FINE = {
    0.850: (0.0728, 0.0358),
    0.875: (0.0743, 0.0401),
    0.900: (0.0774, 0.0462),
    0.925: (0.0799, 0.0529),
    0.950: (0.0800, 0.0565),
    0.975: (0.0727, 0.0510),
    1.000: (0.0628, 0.0419),
    1.025: (0.0590, 0.0381),
    1.050: (0.0579, 0.0375),
    1.075: (0.0579, 0.0389),
    1.100: (0.0560, 0.0376),
}

_MACH = sorted({m for m, _ in _ANCHOR})
_ALPHA = [0.0, 1.25, 4.0, 8.0]  # alpha nodes of the M0.8 row


def _bilinear(M, a, field):
    """Bilinear interpolation on the anchor grid. NaN outside.

    Mach slice at a=1.25 uses the fine-mesh campaign (86k cells) when
    available; alpha shape comes from the legacy M0.8 row (only old rows
    have a=0/4/8). shock_x/cpmin legacy-only (not re-extracted in the
    fine campaign).
    """
    idx = {'CL': 0, 'CD': 1, 'shock': 2, 'cpmin': 3}[field]
    Ms = sorted({m for m, _ in _ANCHOR})
    fine = _ANCHOR_FINE
    if field in ('CL', 'CD') and min(fine) <= M <= max(fine):
        fi = 0 if field == 'CL' else 1
        fm = sorted(fine)
        i = max(j for j, m in enumerate(fm) if m <= M)
        i = min(i, len(fm) - 2)
        t = (M - fm[i]) / (fm[i + 1] - fm[i])
        v125 = (1 - t) * fine[fm[i]][fi] + t * fine[fm[i + 1]][fi]
        # alpha shape from the legacy table (a=0/1.25/4/8 at M0.8 row +
        # Mach series at a=1.25): use the old value ratio
        As = [0.0, 1.25, 8.0]
        jd = max(j for j, aa in enumerate(As) if aa <= a)
        jd = min(jd, len(As) - 2)
        s = (a - As[jd]) / (As[jd + 1] - As[jd])
        b0 = _ANCHOR[(0.8, As[jd])][idx]
        b1 = _ANCHOR[(0.8, As[jd + 1])][idx]
        ref = _ANCHOR[(0.8, 1.25)][idx]
        if None in (b0, b1, ref) or abs(ref) < 1e-12:
            return float('nan')
        return float(v125 * ((1 - s) * b0 + s * b1) / ref)
    # legacy route
    As = [0.0, 1.25, 8.0]
    if not (Ms[0] <= M <= Ms[-1] and As[0] <= a <= As[-1]):
        return float('nan')
    def row(Mm, aa):
        return _ANCHOR.get((Mm, aa), (None,) * 4)[idx]
    i = max(i for i, m in enumerate(Ms) if m <= M)
    i = min(i, len(Ms) - 2)
    m0, m1 = Ms[i], Ms[i + 1]
    t = (M - m0) / max(m1 - m0, 1e-12)
    v0, v1 = row(m0, 1.25), row(m1, 1.25)
    if v0 is None or v1 is None:
        return float('nan')
    v_alpha125 = (1 - t) * v0 + t * v1
    j = max(j for j, aa in enumerate(As) if aa <= a)
    j = min(j, len(As) - 2)
    a0, a1 = As[j], As[j + 1]
    s = (a - a0) / max(a1 - a0, 1e-12)
    b0, b1 = row(0.8, a0), row(0.8, a1)
    ref125 = row(0.8, 1.25)
    if None in (b0, b1, ref125) or abs(ref125) < 1e-12:
        return float('nan')
    shape = ((1 - s) * b0 + s * b1) / ref125
    return float(v_alpha125 * shape)


def shock_region_from_anchors(Mach, alpha_deg):
    """Shock-region estimate for the M0.85-1.1 band (anchor-based).

    Returns dict(x_start, x_end, x_shock, M1, cdw, cpmin, source).
    The OF slabs are the ground truth; pocket extent is derived from
    the local-Mach inversion of the anchor upper Cp (where available)
    or bracketed from the anchor shock location when not.
    """
    Mach = float(Mach); alpha_deg = float(alpha_deg)
    Ms = sorted({m for m, _ in _ANCHOR})
    m_lo = min(Ms[0], min(_ANCHOR_FINE))
    m_hi = max(Ms[-1], max(_ANCHOR_FINE))
    if not (m_lo <= Mach <= m_hi and 0.0 <= alpha_deg <= 8.0):
        return {'x_start': None, 'x_end': None, 'x_shock': None,
                'M1': None, 'cdw': None, 'cpmin': None,
                'source': 'anchor (out of range)'}
    cl = _bilinear(Mach, alpha_deg, 'CL')
    cd = _bilinear(Mach, alpha_deg, 'CD')
    sh = _bilinear(Mach, alpha_deg, 'shock')
    cpm = _bilinear(Mach, alpha_deg, 'cpmin')
    # wave drag = CD - subsonic floor (~0.044 at alpha=1.25)
    floor = 0.044 if alpha_deg <= 2 else 0.040
    cdw = max(cd - floor, 0.0)
    # pocket brackets the shock: start ~ x_shock - 0.25 (typical extent
    # from OF M0.9: 0.11->0.56 with shock 0.40); end ~ shock + 0.10
    if sh is None:
        xstart = xend = None
    else:
        xstart = max(sh - 0.28, 0.05)
        xend = min(sh + 0.22, 0.98)
    return {'x_start': xstart, 'x_end': xend, 'x_shock': sh,
            'M1': None, 'cdw': cdw, 'cpmin': cpm,
            'source': 'anchor'}


def predict_full(Mach, alpha_deg):
    """Section prediction + shock region for M0.65-1.1 (surrogate)."""
    Mach = float(Mach); alpha_deg = float(alpha_deg)
    base = predict(Mach, alpha_deg)
    if base['CL'] is None or not np.isfinite(base['CL']):
        return dict(base, x_start=None, x_end=None, cdw=None)
    reg = shock_region_from_anchors(Mach, alpha_deg)
    base.update(reg)
    return base


def predict(Mach, alpha_deg):
    """Transonic section prediction. Returns dict(CL, CD, shock_x,
    cpmin_up, source). source='anchor' (M0.65-1.1) or 'KT' (<0.65)."""
    Mach = float(Mach)
    alpha_deg = float(alpha_deg)
    if Mach < 0.65:
        return {'CL': float('nan'), 'CD': float('nan'),
                'shock_x': None, 'cpmin_up': None,
                'source': 'KT (use panel+karman_tsien)'}
    Ms = sorted({m for m, _ in _ANCHOR})
    m_lo = min(Ms[0], min(_ANCHOR_FINE))
    m_hi = max(Ms[-1], max(_ANCHOR_FINE))
    if not (m_lo <= Mach <= m_hi and 0.0 <= alpha_deg <= 8.0):
        return {'CL': float('nan'), 'CD': float('nan'),
                'shock_x': None, 'cpmin_up': None,
                'source': 'anchor (out of range)'}
    return {'CL': _bilinear(Mach, alpha_deg, 'CL'),
            'CD': _bilinear(Mach, alpha_deg, 'CD'),
            'shock_x': _bilinear(Mach, alpha_deg, 'shock'),
            'cpmin_up': _bilinear(Mach, alpha_deg, 'cpmin'),
            'source': 'anchor'}


if __name__ == '__main__':
    for M in [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1]:
        p = predict(M, 1.25)
        print('M=%.2f CL=%.3f CD=%.4f shock=%s' % (
            M, p['CL'], p['CD'], p['shock_x']))
    for a in [0, 2, 4, 6, 8]:
        p = predict(0.8, a)
        print('a=%d CL=%.3f CD=%.4f' % (a, p['CL'], p['CD']))


def cp_critical(Mach, gamma=1.4):
    """Critical (sonic) pressure coefficient at freestream Mach."""
    M2 = Mach ** 2
    t = (1 + 0.5 * (gamma - 1) * M2) / (1 + 0.5 * (gamma - 1))
    return 2.0 / (gamma * M2) * (t ** (gamma / (gamma - 1)) - 1.0)


def karman_tsien_cp(cp0, Mach):
    """KT compressibility correction (scalar)."""
    b = max(1 - Mach ** 2, 1e-9)
    return cp0 / (b + Mach ** 2 / (1 + b) * cp0 / 2.0)


def m_crit(cp0min, gamma=1.4):
    """Critical Mach: suction peak first reaches sonic.

    PG-amplified suction vs exact isentropic Cp* intersection
    (bisection). PG matches NACA0012 Mcrit data (~0.73 at low alpha);
    KT-based estimate would sit ~0.07 lower (KT over-deepens --
    conservative, noted but not used).
    Below M_crit the flow is shock-free and KT/PG are defensible;
    above it pockets/shocks form (need TSD/OF/surrogate).
    """
    import numpy as np
    lo, hi = 0.05, 0.99
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        pg = cp0min / max(np.sqrt(1 - mid ** 2), 1e-9)
        if pg < cp_critical(mid, gamma):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def kt_status(Mach, cp0min):
    """Validity verdict for KT/PG at (Mach, suction peak)."""
    mc = m_crit(cp0min)
    if Mach < mc - 0.03:
        return 'valid (subcritical, margin %.2f)' % (mc - Mach)
    if Mach < mc + 0.02:
        return 'marginal (near M_crit=%.2f)' % mc
    return 'invalid (supercritical, M_crit=%.2f)' % mc
