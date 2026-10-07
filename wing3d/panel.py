"""Full-3D Morino Dirichlet panel solver (Hess & Smith family).

Constant source + constant doublet on arbitrary planar polygons, frozen
wake with Morino-Kutta condition. Conventions (self-consistent):
  sigma_j = -Vinf . n_j (known)
  doublet potential  Phi_d = -mu * Omega / (4pi), Omega = signed solid angle
                      (+ on the outer side)
  Dirichlet: interior potential = freestream potential
      sum_j C_ij mu_j + wake = -sum_j B_ij sigma_j,  C = -Omega/4pi
  Kutta per strip: mu_wake = mu_upper_TE - mu_lower_TE
Velocity of a constant doublet panel = vortex loop of strength mu around
its edges (Katz & Plotkin); source velocity by quadrature, self += n/2.
"""
import numpy as np


def _solid_angle_tri(p, a, b, c):
    """Signed solid angle of triangle abc seen from p (van Oosterom-Strackee)."""
    a1 = a - p
    b1 = b - p
    c1 = c - p
    la = np.linalg.norm(a1, axis=-1)
    lb = np.linalg.norm(b1, axis=-1)
    lc = np.linalg.norm(c1, axis=-1)
    num = np.einsum('...i,...i', a1, np.cross(b1, c1))
    den = (la * lb * lc + np.einsum('...i,...i', a1, b1) * lc
           + np.einsum('...i,...i', b1, c1) * la
           + np.einsum('...i,...i', c1, a1) * lb)
    return 2.0 * np.arctan2(num, den)


def solid_angle_poly(p, verts):
    """Signed solid angle of polygon (fan from vertex 0). p: (3,) or (M,3)."""
    p = np.asarray(p, dtype=float)
    v = np.asarray(verts, dtype=float)
    single = (p.ndim == 1)
    p = np.atleast_2d(p)
    om = np.zeros(len(p))
    a = v[0]
    for k in range(1, len(v) - 1):
        om += _solid_angle_tri(p, np.broadcast_to(a, p.shape),
                               np.broadcast_to(v[k], p.shape),
                               np.broadcast_to(v[k + 1], p.shape))
    return om[0] if single else om


def _fan_quads(verts):
    """Quadrature points+weights (area-weighted) over polygon via 4x subdiv fan."""
    v = np.asarray(verts, dtype=float)
    pts, wts = [], []
    a = v[0]
    for k in range(1, len(v) - 1):
        b, c = v[k], v[k + 1]
        m_ab = 0.5 * (a + b)
        m_bc = 0.5 * (b + c)
        m_ca = 0.5 * (c + a)
        g = (a + b + c) / 3.0
        for t in ([a, m_ab, m_ca], [m_ab, b, m_bc],
                  [m_ca, m_bc, c], [m_ab, m_bc, m_ca]):
            t = np.array(t)
            e1 = t[1] - t[0]
            e2 = t[2] - t[0]
            ar = 0.5 * np.linalg.norm(np.cross(e1, e2))
            pts.append(t.mean(axis=0))
            wts.append(ar)
    return np.array(pts), np.array(wts)


def source_potential(p, verts, area):
    """(1/4pi) int 1/r dS. Monopole far-field, quadrature near-field."""
    p = np.asarray(p, dtype=float)
    v = np.asarray(verts, dtype=float)
    c = v.mean(axis=0)
    r = np.linalg.norm(p - c)
    size = np.sqrt(area)
    if r > 4 * size:
        return area / (4 * np.pi * r)
    q, w = _fan_quads(v)
    d = np.linalg.norm(q - p, axis=1)
    d = np.maximum(d, 1e-12 * (1 + size))
    return float(np.sum(w / d) / (4 * np.pi))


def source_velocity(p, verts, area, self_term=False, normal=None):
    """-(1/4pi) int r_vec/r^3 dS; self_term adds exterior +n/2 jump."""
    p = np.asarray(p, dtype=float)
    v = np.asarray(verts, dtype=float)
    out = np.zeros(3)
    if not self_term:
        c = v.mean(axis=0)
        r = np.linalg.norm(p - c)
        size = np.sqrt(area)
        if r > 4 * size:
            rv = p - c
            return -area * rv / (4 * np.pi * r**3)
        q, w = _fan_quads(v)
        rv = q - p
        d = np.linalg.norm(rv, axis=1)
        d = np.maximum(d, 1e-12 * (1 + size))
        out = -np.sum((w / d**3)[:, None] * rv, axis=0) / (4 * np.pi)
    else:
        q, w = _fan_quads(v)
        rv = q - p
        d = np.linalg.norm(rv, axis=1)
        d = np.maximum(d, 1e-12 * (1 + np.sqrt(area)))
        out = -np.sum((w / d**3)[:, None] * rv, axis=0) / (4 * np.pi)
        out += 0.5 * np.asarray(normal, dtype=float)
    return out


def biot_savart(p, e0, e1, gamma=1.0, cutoff=1e-10):
    """Velocity at p from straight vortex segment e0->e1, strength gamma."""
    r0 = e1 - e0
    r1 = p - e0
    r2 = p - e1
    cr = np.cross(r1, r2)
    denom = np.dot(cr, cr) + cutoff
    d1 = max(np.linalg.norm(r1), 1e-14)
    d2 = max(np.linalg.norm(r2), 1e-14)
    cos_term = np.dot(r0, r1 / d1 - r2 / d2)
    return gamma / (4 * np.pi) * cr / denom * cos_term


def loop_velocity(p, verts, gamma=1.0):
    """Velocity of closed vortex loop (CCW seen from tip of normal)."""
    v = np.asarray(verts, dtype=float)
    out = np.zeros(3)
    m = len(v)
    for k in range(m):
        out += biot_savart(p, v[k], v[(k + 1) % m], gamma)
    return out
