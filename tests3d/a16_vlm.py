"""A16: non-planar vortex-lattice reference (pure Python, no container).

Horseshoe VLM on an arbitrary spanwise station list, so the surface may be
planar, dihedraled, tapered, swept, twisted or cranked.  This is the
"defensible non-planar reference" for the planarity-error table; it is a
standard small-disturbance method (exact for the linear lifting problem it
solves), NOT a validated transonic reference.

Formulation
    panel k between stations k,k+1: corners A=LE(k), B=LE(k+1),
    C=TE(k+1), D=TE(k);  normal n = (B-A)x(D-A) normalised, +z sense.
    bound vortex at the quarter-chord line, collocation at mid span /
    3/4 chord, horseshoe legs trailed to x = XFAR (200 chords: the far part
    of the legs contributes O(1/r) at the wing, checked converged).
    A.Gamma = -(Vinf . n)  ->  Gamma per panel.
    lift: f = rho*Gamma*Vinf x dl  =>  f_z = rho*V*Gamma*dy*cos(dihedral),
    i.e. CL = 2*sum(Gamma*dy*cosgam)/S  (the classic cos(dihedral) factor
    falls out of the cross product, it is not put in by hand).

  python a16_vlm.py            validate + print the reference table
"""
import sys

import numpy as np

XFAR = 200.0


def _seg(p, p1, p2):
    """Velocity induced at p by a straight vortex segment p1->p2 of unit
    circulation (Biot-Savart, vectorised over p)."""
    r1 = p - p1
    r2 = p - p2
    a = np.linalg.norm(r1, axis=-1)
    b = np.linalg.norm(r2, axis=-1)
    cross = np.cross(r1, r2)
    cc = np.einsum('...i,...i->...', cross, cross)
    den = a * b * (a * b + np.einsum('...i,...i->...', r1, r2))
    with np.errstate(divide='ignore', invalid='ignore'):
        u = ((a[:, None] + b[:, None]) * cross) / (4 * np.pi * den[:, None])
    bad = (~np.isfinite(u)).any(axis=-1) | (den < 1e-14)
    if np.any(bad):
        # straight-through limit: u -> 0 on the segment, but a collocation
        # point must never coincide with a filament; flag loudly instead
        raise FloatingPointError('collocation point on a vortex segment '
                                 '(%d points)' % int(bad.sum()))
    return u


def stations_rect(span=1.6, chord=1.0, ns=24, xle=0.0, zfun=None,
                  taper=1.0, sweep_deg=0.0, twist_deg=0.0,
                  breakpoints=None):
    """Station list for a (tapered/swept/twisted/non-planar) wing.

    y runs -span/2..+span/2; xle and chord follow taper/sweep; zfun(y)
    gives the surface height (dihedral etc.); twist is linear in eta and
    positive = nose up.
    """
    y = np.linspace(-0.5 * span, 0.5 * span, ns + 1)
    eta = 2.0 * np.abs(y) / span
    c = chord * (1.0 - (1.0 - taper) * eta)
    xl = xle + np.tan(np.radians(sweep_deg)) * np.abs(y)
    z = np.zeros_like(y) if zfun is None else np.asarray(zfun(y), float)
    tw = np.radians(twist_deg) * (1.0 - eta)   # root twist -> 0 at tip
    return {'y': y, 'c': c, 'xle': xl, 'z': z, 'tw': tw}


def vlm(st, alpha_deg=1.25, M=0.0, pg=True):
    """Solve the lattice. Returns dict(CL, Gamma, y, cp...)"""
    y = np.asarray(st['y'], float)
    c = np.asarray(st['c'], float)
    xl = np.asarray(st['xle'], float)
    z = np.asarray(st['z'], float)
    tw = np.asarray(st['tw'], float)
    ns = len(y) - 1
    LE = np.stack([xl, y, z], axis=-1)
    TE = np.stack([xl + c * np.cos(tw), y, z - c * np.sin(tw)], axis=-1)
    A = LE[:-1]
    B = LE[1:]
    C = TE[1:]
    D = TE[:-1]
    # panel normal (+z sense)
    n = np.cross(B - A, D - A)
    n *= np.sign(n[:, 2])[:, None]
    n /= np.linalg.norm(n, axis=-1)[:, None]
    cp = A + 0.5 * (B - A) + 0.75 * (D - A)           # mid span, 3/4 chord
    P1 = A + 0.25 * (D - A)                            # quarter chord, y_k
    P2 = B + 0.25 * (C - B)                            # quarter chord, y_k+1
    far1 = np.stack([np.full(ns, XFAR), P1[:, 1], P1[:, 2]], axis=-1)
    far2 = np.stack([np.full(ns, XFAR), P2[:, 1], P2[:, 2]], axis=-1)
    Aman = np.empty((ns, ns))
    for k in range(ns):
        u = (_seg(cp, far1[k], P1[k]) + _seg(cp, P1[k], P2[k])
             + _seg(cp, P2[k], far2[k]))
        Aman[:, k] = np.einsum('ij,ij->i', u, n)
    al = np.radians(alpha_deg)
    V = np.array([np.cos(al), 0.0, np.sin(al)])
    beta = np.sqrt(max(1.0 - M ** 2, 1e-6)) if pg else 1.0
    rhs = -np.einsum('ij,j->i', n, V) / beta     # Prandtl-Glauert on the
    G = np.linalg.solve(Aman, rhs)               # normal-flow equation
    dy = np.diff(y)
    sg = np.stack([np.zeros(ns), dy, np.diff(z)], axis=-1)
    sg /= np.linalg.norm(sg, axis=-1)[:, None]
    S = float(np.sum(0.5 * (c[:-1] + c[1:]) * np.abs(dy)))
    CL = 2.0 * float(np.sum(G * dy * sg[:, 1])) / S
    return {'CL': CL, 'Gamma': G, 'y': y, 'c': c, 'z': z, 'S': S,
            'cl_sec': 2.0 * G / np.maximum(c[:-1], 1e-12),
            'eta': np.abs(y) / (0.5 * (y[-1] - y[0]))}


def _ll(alpha_deg, span, chord=1.0, M=0.0):
    cl2 = 2 * np.pi * np.radians(alpha_deg)
    beta = max(np.sqrt(max(1 - M ** 2, 1e-6)), 0.3)
    AR = span / chord
    return cl2 / beta * AR / (AR + 2.0)


def validate():
    print('--- VLM validation (planar, no container needed) ---')
    for ns in (12, 24, 48, 96):
        for AR in (6.0, 1.6):
            st = stations_rect(span=AR, ns=ns)
            r = vlm(st, 1.25, 0.0)
            ll = _ll(1.25, AR)
            print('  ns=%3d AR=%.1f  CL=%.5f   LL(AR/(AR+2))=%.5f  '
                  'ratio=%.4f' % (ns, AR, r['CL'], ll, r['CL'] / ll))
    print('\n--- dihedral: is dCL = cos(Gamma)? ---')
    st0 = stations_rect(span=1.6, ns=48)
    c0 = vlm(st0, 1.25, 0.0)['CL']
    for g in (0.0, 2.0, 5.0, 10.0, 20.0):
        st = stations_rect(span=1.6, ns=48,
                           zfun=lambda yy, g=g: np.abs(yy) * np.tan(
                               np.radians(g)))
        cl = vlm(st, 1.25, 0.0)['CL']
        print('  dih %5.1f deg  CL=%.5f  dCL=%+.4f%%  cos=%.5f  '
              'dCL/cos-1=%+.4f%%' % (
                  g, cl, 100 * (cl - c0) / c0, np.cos(np.radians(g)),
                  100 * ((cl - c0) / c0 + 1 - np.cos(np.radians(g)))))
    print('\n--- spanwise height variation (non-dihedral) ---')
    for h in (0.05, 0.10, 0.20):
        st = stations_rect(
            span=1.6, ns=48,
            zfun=lambda yy, h=h: h * np.cos(np.pi * yy / 1.6))
        cl = vlm(st, 1.25, 0.0)['CL']
        slope2 = (np.pi * h / 1.6) ** 2
        print('  arch h=%.2fc  CL=%.5f  dCL=%+.4f%%   bound '
              '0.5*(dh/dy)^2=%+.4f%%' % (
                  h, cl, 100 * (cl - c0) / c0, 100 * 0.5 * slope2))
    for h in (0.05, 0.10):
        st = stations_rect(
            span=1.6, ns=48,
            zfun=lambda yy, h=h: h * (1.0 - (2 * yy / 1.6) ** 2) ** 2)
        cl = vlm(st, 1.25, 0.0)['CL']
        print('  bump h=%.2fc  CL=%.5f  dCL=%+.4f%%' % (
            h, cl, 100 * (cl - c0) / c0))


if __name__ == '__main__':
    validate()
