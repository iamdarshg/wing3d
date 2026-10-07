"""Minimal vortex-ring VLM (Katz & Plotkin) as independent lifting arbiter."""
import numpy as np


def _vort_seg(p, e0, e1, gamma=1.0):
    r0 = e1 - e0
    r1 = p - e0
    r2 = p - e1
    cr = np.cross(r1, r2)
    den = np.dot(cr, cr) + 1e-14
    d1 = max(np.linalg.norm(r1), 1e-14)
    d2 = max(np.linalg.norm(r2), 1e-14)
    return gamma / (4 * np.pi) * cr / den * np.dot(r0, r1 / d1 - r2 / d2)


def vlm_cl(span=6.0, chord=1.0, alpha_deg=4.0, n_span=24, n_chord=6):
    """Flat rectangular wing VLM. Returns CL (V=1, rho=1)."""
    a = np.radians(alpha_deg)
    V = np.array([np.cos(a), 0.0, np.sin(a)])
    ys = np.linspace(-span / 2, span / 2, n_span + 1)
    xs = np.linspace(0, chord, n_chord + 1)
    rings = []
    for j in range(n_span):
        for i in range(n_chord):
            x0, x1 = xs[i], xs[i + 1]
            y0, y1 = ys[j], ys[j + 1]
            # bound at x0+dx/4, control at x0+3dx/4 (1/4-3/4 rule)
            rings.append(dict(xb0=x0, xb1=x1, y0=y0, y1=y1,
                              xc=(x0 + 3 * (x1 - x0) / 4), yc=(y0 + y1) / 2))
    n = len(rings)
    A = np.zeros((n, n))
    rhs = np.zeros(n)
    TE = chord + 5 * chord
    for r, R in enumerate(rings):
        pc = np.array([R['xc'], R['yc'], 0.0])
        rhs[r] = -(V @ np.array([0., 0., 1.]))
        for c, C in enumerate(rings):
            # horseshoe: bound segment + two trailing legs to downstream inf
            b0 = np.array([C['xb0'], C['y0'], 0.0])
            b1 = np.array([C['xb0'], C['y1'], 0.0])
            t0 = np.array([TE, C['y0'], 0.0])
            t1 = np.array([TE, C['y1'], 0.0])
            v = (_vort_seg(pc, b0, b1) + _vort_seg(pc, b1, t1)
                 + _vort_seg(pc, t0, b0))
            A[r, c] = v[2]  # normal (z) velocity
    gam = np.linalg.solve(A, rhs)
    # forces via Kutta-Joukowski on bound segments
    L = 0.0
    for g, C in zip(gam, rings):
        b0 = np.array([C['xb0'], C['y0'], 0.0])
        b1 = np.array([C['xb0'], C['y1'], 0.0])
        dl = b1 - b0
        # local velocity ~ V (plus induced; use V for simplicity? use full)
        L += 1.0 * np.cross(V, g * dl)[2]
    return L / (0.5 * span * chord)
