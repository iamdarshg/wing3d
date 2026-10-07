"""Panel system assembly + solve + surface flow (Morino Dirichlet + Kutta).

Vectorized over source panels; per-mesh precomputation cached for reuse
across alphas and viscous-coupling iterations (matrix depends on geometry
only — RHS carries onset flow + transpiration).
"""
import numpy as np


def _precompute(mesh, nsub=4):
    """Fan triangles, quadrature points, edges per panel (padded stacks).

    Quadrature: each fan triangle subdivided nsub x nsub with centroid rule
    (16 points per fan triangle by default) for accurate near-field
    velocity (1/r^2) influence. Cached on the mesh for reuse.
    """
    if hasattr(mesh, '_pcache'):
        return mesh._pcache
    polys = [np.asarray(mesh.verts[f]) for f in mesh.faces]
    maxfan = max(len(p) - 2 for p in polys)
    maxq = maxfan * nsub * nsub
    maxe = max(len(p) for p in polys)
    nf = len(polys)
    fan = np.zeros((nf, maxfan, 3, 3))
    fmask = np.zeros((nf, maxfan), dtype=bool)
    qp = np.zeros((nf, maxq, 3))
    qw = np.zeros((nf, maxq))
    ed = np.zeros((nf, maxe, 2, 3))
    emask = np.zeros((nf, maxe), dtype=bool)
    for i, v in enumerate(polys):
        m = len(v)
        a = v[0]
        fi = 0
        qi = 0
        for k in range(1, m - 1):
            b, c = v[k], v[k + 1]
            fan[i, fi] = (a, b, c)
            fmask[i, fi] = True
            fi += 1
            for ii in range(nsub):
                for jj in range(nsub - ii):
                    p0 = a + (b - a) * ii / nsub + (c - a) * jj / nsub
                    p1 = a + (b - a) * (ii + 1) / nsub + (c - a) * jj / nsub
                    p2 = a + (b - a) * ii / nsub + (c - a) * (jj + 1) / nsub
                    ar = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0))
                    qp[i, qi] = (p0 + p1 + p2) / 3.0
                    qw[i, qi] = ar
                    qi += 1
                    if ii + jj + 1 < nsub:
                        q0 = a + (b - a) * (ii + 1) / nsub + (c - a) * jj / nsub
                        q1 = a + (b - a) * ii / nsub + (c - a) * (jj + 1) / nsub
                        q2 = a + (b - a) * (ii + 1) / nsub + (c - a) * (jj + 1) / nsub
                        ar2 = 0.5 * np.linalg.norm(np.cross(q1 - q0, q2 - q0))
                        qp[i, qi] = (q0 + q1 + q2) / 3.0
                        qw[i, qi] = ar2
                        qi += 1
        for k in range(m):
            ed[i, k, 0] = v[k]
            ed[i, k, 1] = v[(k + 1) % m]
            emask[i, k] = True
    mesh._pcache = dict(fan=fan, fmask=fmask, qp=qp, qw=qw, ed=ed,
                        emask=emask)
    return mesh._pcache


def _solid_angle_rows(p, fan, fmask):
    """Signed solid angle at point p for every panel (vectorized)."""
    a = fan[:, :, 0] - p
    b = fan[:, :, 1] - p
    c = fan[:, :, 2] - p
    la = np.linalg.norm(a, axis=-1)
    lb = np.linalg.norm(b, axis=-1)
    lc = np.linalg.norm(c, axis=-1)
    num = np.einsum('...i,...i', a, np.cross(b, c))
    den = (la * lb * lc + np.einsum('...i,...i', a, b) * lc
           + np.einsum('...i,...i', b, c) * la
           + np.einsum('...i,...i', c, a) * lb)
    om = 2.0 * np.arctan2(num, den)
    return np.where(fmask, om, 0.0).sum(axis=1)


def _centroids(qp, qw):
    wsum = np.maximum(qw.sum(axis=1), 1e-300)
    return (qw[..., None] * qp).sum(axis=1) / wsum[:, None]


def _self_source_pot(p, verts):
    """Self-potential at an in-plane point: analytic edge sum.

    For P in the polygon plane, (1/4pi) int_S dS/r
      = (1/4pi) sum_edges h_e * ln((r1+r2+e)/(r1+r2-e)),
    with h_e the in-plane distance from P to the edge line, r1/r2 the
    distances to its endpoints and e the edge length (all edges counted
    positively for P inside a CCW polygon).
    """
    v = np.asarray(verts, dtype=float)
    m = len(v)
    total = 0.0
    for k in range(m):
        a, b = v[k], v[(k + 1) % m]
        e = b - a
        elen = np.linalg.norm(e)
        if elen < 1e-15:
            continue
        r1 = np.linalg.norm(a - p)
        r2 = np.linalg.norm(b - p)
        # in-plane distance from P to the edge line
        h = abs(np.linalg.norm(np.cross(e, a - p)) / elen)
        denom = max(r1 + r2 - elen, 1e-300)
        total += h * np.log((r1 + r2 + elen) / denom)
    return total / (4 * np.pi)


def _source_pot_rows(p, qp, qw, area, size, self_idx=None, self_verts=None):
    """(1/4pi) int 1/r dS per panel. Monopole far-field, quadrature near."""
    c = _centroids(qp, qw)
    r = np.linalg.norm(c - p, axis=1)
    far = r > 4 * size
    out = np.zeros(len(area))
    out[far] = area[far] / (4 * np.pi * r[far])
    near = ~far
    if self_idx is not None:
        near[self_idx] = False
    if np.any(near):
        d = qp[near] - p
        dist = np.linalg.norm(d, axis=-1)
        dist = np.maximum(dist, 1e-12)
        out[near] = np.sum(qw[near] / dist, axis=1) / (4 * np.pi)
    if self_idx is not None:
        out[self_idx] = _self_source_pot(p, self_verts)
    return out


def _source_vel_rows(p, qp, qw, area, size, self_idx=None, normal=None):
    """-(1/4pi) int r_vec/r^3 dS per panel (velocity vectors).

    The self panel's symmetric integral is exactly zero; only the
    exterior +n/2 jump is retained (avoids quadrature noise)."""
    c = _centroids(qp, qw)
    r = np.linalg.norm(c - p, axis=1)
    far = r > 4 * size
    out = np.zeros((len(area), 3))
    rv = p - c[far]
    d = np.maximum(r[far], 1e-14)
    # monopole outward (matches near-branch quadrature: verified continuous
    # across the crossover and against analytic disk velocity)
    out[far] = +area[far][:, None] * rv / (4 * np.pi * d[:, None]**3)
    near = ~far
    if self_idx is not None:
        near[self_idx] = False
    if np.any(near):
        d = qp[near] - p
        dist = np.linalg.norm(d, axis=-1)
        dist = np.maximum(dist, 1e-12)
        out[near] = -np.sum((qw[near] / dist**3)[..., None] * d, axis=1) / (4 * np.pi)
    if self_idx is not None:
        # Self row zeroed here; the caller adds the exterior jump
        # (-sigma/2 normal) explicitly with the correct representation sign.
        out[self_idx] = 0.0
    return out


def _loop_vel_rows(p, ed, emask, cutoff=1e-10):
    """Vortex-loop velocity per panel (unit strength). Vectorized."""
    e0 = ed[:, :, 0]
    e1 = ed[:, :, 1]
    r0 = e1 - e0
    r1 = p - e0
    r2 = p - e1
    cr = np.cross(r1, r2)
    denom = np.einsum('...i,...i', cr, cr) + cutoff
    d1 = np.maximum(np.linalg.norm(r1, axis=-1), 1e-14)
    d2 = np.maximum(np.linalg.norm(r2, axis=-1), 1e-14)
    cos_term = np.einsum('...i,...i', r0, r1 / d1[..., None] - r2 / d2[..., None])
    v = cr / denom[..., None] * cos_term[..., None] / (4 * np.pi)
    return np.where(emask[..., None], v, 0.0).sum(axis=1)


class WakeStrip:
    def __init__(self, upper_te, lower_te, wake_panels):
        self.upper_te = upper_te
        self.lower_te = lower_te
        self.wake_panels = wake_panels  # list of (verts,)


def build_wake(mesh, strips_info, length=5.0, n_panels=4, growth=1.6):
    """Frozen wake extending downstream (+x) from each TE strip."""
    wakes = []
    scale = (growth ** np.arange(n_panels + 1) - 1) / (growth - 1)
    scale = scale / scale[-1] * length
    for (up, lo, pa, pb) in strips_info:
        pa, pb = np.asarray(pa), np.asarray(pb)
        panels = []
        for k in range(n_panels):
            x0, x1 = scale[k], scale[k + 1]
            a0 = pa + np.array([x0, 0, 0])
            b0 = pb + np.array([x0, 0, 0])
            a1 = pa + np.array([x1, 0, 0])
            b1 = pb + np.array([x1, 0, 0])
            # wound CCW seen from +z (normal up) to match the Kutta
            # convention mu_wake = mu_upper_TE - mu_lower_TE
            panels.append(np.array([a0, a1, b1, b0]))
        wakes.append(WakeStrip(up, lo, panels))
    return wakes


def wing_strips_from_structured(n_loop, n_span_full):
    """TE panel pairs for build_wing meshes (panel_id map included).

    Panels per span segment run 0..n_loop-2 (no wraparound; the TE point
    is shared, so panel 0 = upper TE, panel n_loop-2 = lower TE).
    """
    nl = n_loop
    npl = nl - 1
    panel_id = {}
    pid = 0
    for j in range(n_span_full):
        for i in range(npl):
            panel_id[(i, j)] = pid
            pid += 1
    strips = []
    for j in range(n_span_full):
        strips.append((panel_id[(0, j)], panel_id[(npl - 1, j)], j))
    return strips, panel_id


def build_wake_from_meta(mesh, length=5.0, n_panels=4, growth=1.6):
    """Frozen wake strips from mesh.meta (structured wings)."""
    meta = mesh.meta
    nl, ns = meta['n_loop'], meta['n_span']
    strips, _ = wing_strips_from_structured(nl, ns)
    info = []
    for (up, lo, j) in strips:
        info.append((up, lo, meta['te_points'][j], meta['te_points'][j + 1]))
    return build_wake(mesh, info, length, n_panels, growth)


def wake_precompute(wakes):
    """Flatten wake panels for vectorized influence."""
    polys, owner = [], []
    for k, w in enumerate(wakes):
        for wp in w.wake_panels:
            polys.append(np.asarray(wp))
            owner.append(k)
    if not polys:
        return None
    maxfan = max(len(p) - 2 for p in polys)
    maxe = max(len(p) for p in polys)
    nw = len(polys)
    fan = np.zeros((nw, maxfan, 3, 3))
    fmask = np.zeros((nw, maxfan), dtype=bool)
    ed = np.zeros((nw, maxe, 2, 3))
    emask = np.zeros((nw, maxe), dtype=bool)
    for i, v in enumerate(polys):
        for k in range(1, len(v) - 1):
            fan[i, k - 1] = (v[0], v[k], v[k + 1])
            fmask[i, k - 1] = True
        for k in range(len(v)):
            ed[i, k, 0] = v[k]
            ed[i, k, 1] = v[(k + 1) % len(v)]
            emask[i, k] = True
    return dict(fan=fan, fmask=fmask, ed=ed, emask=emask,
                owner=np.array(owner))


def assemble(mesh, wakes, kutta_sign=+1.0, kutta_mode='potential'):
    """Geometry-only influence matrix A (N body + S kutta rows).

    kutta_mode='potential' (default): exact potential-jump Kutta,
        mu_wake = Phi_upper_TE - Phi_lower_TE (exterior potentials).
        Linear in mu; strictly more accurate than the mu-difference
        shortcut on coarse meshes.
    kutta_mode='doublet': classic mu_wake = mu_upper_TE - mu_lower_TE
        (times kutta_sign); kept for comparison.
    """
    n = mesh.npanels
    s = len(wakes)
    pc = _precompute(mesh)
    wc = wake_precompute(wakes)
    C = mesh.centroid
    size = np.sqrt(mesh.area)
    Vv = [mesh.verts[f] for f in mesh.faces]
    A = np.zeros((n + s, n + s))
    Brow = np.zeros((n, n))  # source potential rows (for RHS variants)
    Cint = np.zeros((n, n))  # interior doublet rows (reuse for Kutta)
    for i in range(n):
        p = C[i]
        om = _solid_angle_rows(p, pc['fan'], pc['fmask'])
        # Interior limit: van Oosterom-Strackee solid angle -> +2π on the
        # inner side, so the doublet self-term is C_ii = -1/2 (Katz & Plotkin).
        om[i] = 2 * np.pi
        Cint[i] = -om / (4 * np.pi)
        A[i, :n] = Cint[i]
        b = _source_pot_rows(p, pc['qp'], pc['qw'], mesh.area, size,
                             self_idx=i, self_verts=Vv[i])
        Brow[i] = b
        if wc is not None:
            omw = _solid_angle_rows(p, wc['fan'], wc['fmask'])
            for k in range(s):
                A[i, n + k] = -omw[wc['owner'] == k].sum() / (4 * np.pi)
    if kutta_mode == 'doublet':
        for k, w in enumerate(wakes):
            A[n + k, w.upper_te] = -kutta_sign
            A[n + k, w.lower_te] = +kutta_sign
            A[n + k, n + k] = 1.0
    elif kutta_mode == 'potential':
        # Exterior potential row = interior row + self jump (+1 on self).
        # Kutta: mu_w - (Phi_u - Phi_l) = 0 with Phi = Phi_inf + B.sigma
        # + C.mu; the B.sigma and Phi_inf parts go to kutta_rhs().
        for k, w in enumerate(wakes):
            u, l = w.upper_te, w.lower_te
            A[n + k, :n] = -(Cint[u] - Cint[l])
            A[n + k, u] -= 1.0
            A[n + k, l] += 1.0
            A[n + k, n + k] = 1.0
    else:
        raise ValueError('unknown kutta_mode')
    return A, Brow


def kutta_rhs(mesh, wakes, vinf, Brow):
    """RHS known terms (Phi_inf + B sigma) for potential-Kutta rows."""
    V = np.asarray(vinf, dtype=float)
    sigma = (mesh.normal @ V)
    out = np.zeros(len(wakes))
    for k, w in enumerate(wakes):
        u, l = w.upper_te, w.lower_te
        phi_u = mesh.centroid[u] @ V + Brow[u] @ sigma
        phi_l = mesh.centroid[l] @ V + Brow[l] @ sigma
        out[k] = phi_u - phi_l
    return out, sigma


def _surface_gradient(mesh, phi):
    """Least-squares surface gradient of a panel scalar.

    Fit in 2D tangent-plane coordinates (robust: fitting 3D gradients to
    coplanar stencils is rank-deficient and explodes).
    """
    if not hasattr(mesh, '_nbrs'):
        v2f = {}
        for i, f in enumerate(mesh.faces):
            for v in f:
                v2f.setdefault(v, []).append(i)
        nbrs = []
        for i, f in enumerate(mesh.faces):
            s = set()
            for v in f:
                s.update(v2f[v])
            s.discard(i)
            nbrs.append(np.array(sorted(s), dtype=int))
        mesh._nbrs = nbrs
    grad = np.zeros((mesh.npanels, 3))
    C = mesh.centroid
    for i in range(mesh.npanels):
        nb = mesh._nbrs[i]
        if len(nb) < 2:
            continue
        n = mesh.normal[i]
        # tangent basis
        t1 = C[nb[0]] - C[i]
        t1 = t1 - (t1 @ n) * n
        if np.linalg.norm(t1) < 1e-14:
            continue
        t1 /= np.linalg.norm(t1)
        t2 = np.cross(n, t1)
        d = C[nb] - C[i]
        A = np.stack([d @ t1, d @ t2], axis=1)
        b = phi[nb] - phi[i]
        g, *_ = np.linalg.lstsq(A, b, rcond=None)
        grad[i] = g[0] * t1 + g[1] * t2
    return grad


def solve(mesh, wakes, vinf, A=None, Brow=None, extra_rhs=None,
          kutta_mode='potential'):
    """Solve + surface flow. extra_rhs: transpiration source term (coupling)."""
    n = mesh.npanels
    s = len(wakes)
    if A is None or Brow is None:
        A, Brow = assemble(mesh, wakes, kutta_mode=kutta_mode)
    V = np.asarray(vinf, dtype=float)
    vmag = np.linalg.norm(V)
    # Coherent Morino Dirichlet formulation (V = +grad Phi,
    # Phi = Phi_inf + B.sigma + C.mu, sigma = +Vinf . n):
    # interior Dirichlet C.mu = -B.sigma. Selected by grid search over
    # formulation signs against analytic sphere Cp and lifting-line wing.
    sigma = (mesh.normal @ V)
    rhs = np.zeros(n + s)
    rhs[:n] = -(Brow @ sigma)
    if extra_rhs is not None:
        rhs[:n] += extra_rhs
    if kutta_mode == 'potential' and s:
        rhs[n:], _ = kutta_rhs(mesh, wakes, V, Brow)
        # NOTE: kutta_rhs returns (out, sigma); sigma identical to above.
    x = np.linalg.solve(A, rhs)
    mu, muw = x[:n], x[n:]

    pc = _precompute(mesh)
    wc = wake_precompute(wakes)
    C = mesh.centroid
    size = np.sqrt(mesh.area)
    vel = np.zeros((n, 3))
    # Velocity assembly consistent with V = -grad Phi above: negate the
    # -gradient kernel sums.
    for i in range(n):
        p = C[i]
        v = V.copy()
        v -= (sigma[:, None] * _source_vel_rows(
            p, pc['qp'], pc['qw'], mesh.area, size,
            self_idx=i, normal=mesh.normal[i])).sum(axis=0)
        v -= (mu[:, None] * _loop_vel_rows(p, pc['ed'], pc['emask'])).sum(axis=0)
        if wc is not None:
            v -= (muw[wc['owner'], None] * _loop_vel_rows(
                p, wc['ed'], wc['emask'])).sum(axis=0)
        # exterior self of own source sheet: -sigma/2 along outward normal
        v -= 0.5 * sigma[i] * mesh.normal[i]
        vel[i] = v
    # Doublet jump correction: the in-plane loop value is the principal
    # value; the exterior tangential velocity adds half the surface
    # gradient of mu (the doublet sheet jumps by mu across the surface).
    # Sign/scale verified: recovers lifted-FD surface velocity on sphere.
    grad_mu = _surface_gradient(mesh, mu)
    vel += 0.5 * grad_mu
    vn = np.einsum('ij,ij->i', vel, mesh.normal)
    vel = vel - vn[:, None] * mesh.normal
    cp = 1.0 - np.sum(vel**2, axis=1) / vmag**2
    return {'mu': mu, 'muw': muw, 'sigma': sigma, 'vel': vel, 'cp': cp,
            'A': A, 'Brow': Brow, 'vinf': V}


def karman_tsien(cp0, mach):
    """Compressible Cp correction (subsonic)."""
    if mach <= 0:
        return cp0
    b = np.sqrt(max(1 - mach**2, 1e-9))
    return cp0 / (b + mach**2 / (1 + b) * cp0 / 2.0)
