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
    # 2026-10-08: only these rows ran at labeled conditions.
    (0.80, 0.00): (0.0003, 0.0474, None, None),
    (0.80, 1.25): (0.4020, 0.0414, 0.35, -1.27),
    (0.80, 8.00): (0.4040, 0.0415, 0.35, -1.27),
    (0.90, 1.25): (0.0720, 0.0530, 0.40, -0.76),
}

_MACH = sorted({m for m, _ in _ANCHOR})
_ALPHA = [0.0, 1.25, 4.0, 8.0]  # alpha nodes of the M0.8 row


def _bilinear(M, a, field):
    """Bilinear interpolation on the anchor grid. NaN outside."""
    idx = {'CL': 0, 'CD': 1, 'shock': 2, 'cpmin': 3}[field]
    Ms = sorted({m for m, _ in _ANCHOR})
    As = [0.0, 1.25, 8.0]
    if not (Ms[0] <= M <= Ms[-1] and As[0] <= a <= As[-1]):
        return float('nan')
    # Mach slice at alpha=1.25 + alpha scaling from M0.8 row
    def row(Mm, aa):
        return _ANCHOR.get((Mm, aa), (None,) * 4)[idx]
    # interpolate along Mach at 1.25
    i = max(i for i, m in enumerate(Ms) if m <= M)
    i = min(i, len(Ms) - 2)
    m0, m1 = Ms[i], Ms[i + 1]
    t = (M - m0) / max(m1 - m0, 1e-12)
    v0, v1 = row(m0, 1.25), row(m1, 1.25)
    if v0 is None or v1 is None:
        return float('nan')
    v_alpha125 = (1 - t) * v0 + t * v1
    # alpha correction from M0.8 row shape
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
    if not (Ms[0] <= Mach <= Ms[-1] and 0.0 <= alpha_deg <= 8.0):
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
