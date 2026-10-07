"""Forces, moments, Trefftz induced drag. Coefficients use qS with q=0.5*rho*V^2."""
import numpy as np


def pressure_forces(mesh, cp, vinf, rho=1.0, sref=None):
    """Integrate pressure to wind-axis CL/CD/CS + CM. Returns dict."""
    V = np.asarray(vinf, dtype=float)
    vmag = np.linalg.norm(V)
    vhat = V / vmag
    q = 0.5 * rho * vmag ** 2
    F = (-cp[:, None] * mesh.normal * mesh.area[:, None]).sum(axis=0) * q
    # wind axes: drag along vhat, lift perpendicular (in x-z plane), side
    lift_dir = np.array([0.0, 0.0, 1.0])
    lift_dir = lift_dir - (lift_dir @ vhat) * vhat
    lift_dir /= max(np.linalg.norm(lift_dir), 1e-300)
    side_dir = np.cross(vhat, lift_dir)
    side_dir /= max(np.linalg.norm(side_dir), 1e-300)
    drag_dir = np.cross(lift_dir, side_dir)
    if sref is None:
        sref = float(mesh.area.sum()) / 2.0  # wetted/2 ~ planform for wings
    return {'CL': float(F @ lift_dir / (q * sref)),
            'CDp': float(F @ drag_dir / (q * sref)),
            'CS': float(F @ side_dir / (q * sref)),
            'F': F, 'q': q, 'sref': sref,
            'dirs': {'lift': lift_dir, 'drag': drag_dir, 'side': side_dir}}


def spanwise_loading(mesh, cp, vinf, strips_fn, rho=1.0):
    """Sectional cl(y) per spanwise strip. strips_fn maps panel idx -> strip id."""
    V = np.asarray(vinf, dtype=float)
    vmag = np.linalg.norm(V)
    q = 0.5 * rho * vmag ** 2
    vhat = V / vmag
    lift_dir = np.array([0.0, 0.0, 1.0])
    lift_dir = lift_dir - (lift_dir @ vhat) * vhat
    lift_dir /= max(np.linalg.norm(lift_dir), 1e-300)
    out = {}
    for i in range(mesh.npanels):
        k = strips_fn(i)
        d = out.setdefault(k, {'L': 0.0, 'area': 0.0, 'y': []})
        f = -cp[i] * float(mesh.normal[i] @ lift_dir) * mesh.area[i] * q
        d['L'] += f
        d['area'] += mesh.area[i]
        d['y'].append(mesh.centroid[i][1])
    rows = []
    for k, d in sorted(out.items()):
        chord = d['area'] / max(len(d['y']), 1) * 2  # ~chord for 2-sided strip
        rows.append({'strip': k, 'y': float(np.mean(d['y'])),
                     'cl': float(d['L'] / (q * chord * 1.0)) if chord > 0 else 0.0,
                     'chord': float(chord)})
    return rows


def trefftz_cd(cl_span, y_span, vmag=1.0):
    """Induced drag from spanwise cl via Glauert Fourier fit (e <= 1)."""
    y = np.asarray(y_span, dtype=float)
    cl = np.asarray(cl_span, dtype=float)
    b = y.max() - y.min()
    theta = np.arccos(np.clip(2 * (y - y.min()) / max(b, 1e-300) - 1, -1, 1))
def trefftz_cd(cl_span, y_span, chord_span):
    """Induced drag from spanwise cl via Glauert Fourier fit (e <= 1).

    Fits Gamma(y)/(b*V) = sum A_n sin(n*theta) and returns
    CDi = pi*AR*sum(n*A_n^2)/2 with the Oswald efficiency.
    """
    y = np.asarray(y_span, dtype=float)
    cl = np.asarray(cl_span, dtype=float)
    chord = np.asarray(chord_span, dtype=float)
    order = np.argsort(y)
    y, cl, chord = y[order], cl[order], chord[order]
    b = y.max() - y.min()
    theta = np.arccos(np.clip(2 * (y - y.min()) / max(b, 1e-300) - 1, -1, 1))
    gamma = cl * chord / 2.0  # Gamma/V in V=1 units
    S = np.stack([np.sin(n * theta) for n in (1, 3, 5)], axis=1)
    A, *_ = np.linalg.lstsq(S, gamma / max(b, 1e-300), rcond=None)
    area = float(np.trapezoid(chord, y))
    AR = b ** 2 / max(area, 1e-300)
    CDi = float(np.pi * AR * (A[0] ** 2 + 3 * A[1] ** 2 + 5 * A[2] ** 2) / 2)
    CL = float(np.pi * AR * A[0] / 2)
    e = float(CL ** 2 / max(np.pi * AR * CDi, 1e-300)) if CDi > 0 else float('nan')
    return {'CDi': CDi, 'CL': CL, 'e': e, 'A': A}
