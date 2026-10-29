"""Geometry: NACA wing meshes, STL import, decimation. All 3D, no strips."""
import numpy as np

try:
    from scipy.spatial import ConvexHull
    _HAS_QHULL = True
except Exception:
    _HAS_QHULL = False


class Mesh:
    """Polygon surface mesh. verts (Nv,3), faces = list of int lists."""
    def __init__(self, verts, faces):
        self.verts = np.asarray(verts, dtype=float)
        self.faces = [list(map(int, f)) for f in faces]
        self._compute()

    def _compute(self):
        nv = len(self.verts)
        self.centroid = np.zeros((len(self.faces), 3))
        self.normal = np.zeros((len(self.faces), 3))
        self.area = np.zeros(len(self.faces))
        for i, f in enumerate(self.faces):
            p = self.verts[f]
            c = p.mean(axis=0)
            # Newell normal (works for any planar polygon)
            n = np.zeros(3)
            m = len(f)
            for k in range(m):
                a, b = p[k], p[(k + 1) % m]
                n[0] += (a[1] - b[1]) * (a[2] + b[2])
                n[1] += (a[2] - b[2]) * (a[0] + b[0])
                n[2] += (a[0] - b[0]) * (a[1] + b[1])
            area = 0.5 * np.linalg.norm(n)
            self.centroid[i] = c
            self.normal[i] = n / (2 * area) if area > 0 else np.array([0.0, 0.0, 1.0])
            self.area[i] = area
        # orient outward: closed body => sum(c . n * A) > 0
        s = float((self.centroid * self.normal).sum(axis=1) @ self.area)
        if s < 0:
            self.faces = [f[::-1] for f in self.faces]
            self.normal = -self.normal

    @property
    def npanels(self):
        return len(self.faces)

    def write_stl(self, path, binary=True):
        tris = []
        for f in self.faces:
            for k in range(1, len(f) - 1):
                tris.append((f[0], f[k], f[k + 1]))
        tris = np.array(tris, dtype=np.int64)
        v = self.verts
        if binary:
            with open(path, 'wb') as fh:
                fh.write(b'\x00' * 80)
                fh.write(np.uint32(len(tris)))
                for t in tris:
                    n = np.cross(v[t[1]] - v[t[0]], v[t[2]] - v[t[0]])
                    nn = np.linalg.norm(n)
                    n = n / nn if nn > 0 else np.array([0, 0, 1.0])
                    fh.write(np.asarray(n, dtype=np.float32).tobytes())
                    fh.write(np.asarray([v[t[0]], v[t[1]], v[t[2]]],
                                        dtype=np.float32).tobytes())
                    fh.write(np.uint16(0))
        else:
            with open(path, 'w') as fh:
                fh.write('solid wing3d\n')
                for t in tris:
                    fh.write('facet normal 0 0 0\nouter loop\n')
                    for k in t:
                        fh.write(f'vertex {v[k,0]} {v[k,1]} {v[k,2]}\n')
                    fh.write('endloop\nendfacet\n')
                fh.write('endsolid wing3d\n')


def naca4_section(code='0012', n=60, cosine=True):
    """Airfoil loop TE -> upper -> LE -> lower -> TE (no duplicate TE)."""
    code = code.zfill(4)
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0
    if cosine:
        beta = np.linspace(0, np.pi, n + 1)
        x = 0.5 * (1 - np.cos(beta))
    else:
        x = np.linspace(0, 1, n + 1)
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                  + 0.2843 * x**3 - 0.1036 * x**4)  # closed TE
    if m == 0 or p == 0:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    else:
        cond = x < p
        yc = np.where(cond, m / p**2 * (2 * p * x - x**2),
                      m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))
        dyc = np.where(cond, 2 * m / p**2 * (p - x),
                       2 * m / (1 - p)**2 * (p - x))
    th = np.arctan(dyc)
    xu = x - yt * np.sin(th)
    yu = yc + yt * np.cos(th)
    xl = x + yt * np.sin(th)
    yl = yc - yt * np.cos(th)
    upper = np.stack([xu, yu], axis=1)[::-1]   # TE -> LE
    lower = np.stack([xl, yl], axis=1)[1:]     # LE+1 -> TE
    loop = np.vstack([upper, lower])           # closed implicitly
    return loop


def section_loop(xu, yu, xl, yl, n=60):
    """Airfoil loop (TE -> upper -> LE -> lower) from surface coordinates.

    Same convention as naca4_section output: cosine-clustered stations,
    upper TE->LE then lower LE+1->TE. Handles either input ordering.
    """
    xu = np.asarray(xu, float)
    yu = np.asarray(yu, float)
    xl = np.asarray(xl, float)
    yl = np.asarray(yl, float)
    if xu[0] > xu[-1]:
        xu, yu = xu[::-1], yu[::-1]
    if xl[0] > xl[-1]:
        xl, yl = xl[::-1], yl[::-1]
    beta = np.linspace(0.0, np.pi, n + 1)
    xg = 0.5 * (1.0 - np.cos(beta))
    yu_g = np.interp(xg, xu, yu)
    yl_g = np.interp(xg, xl, yl)
    upper = np.stack([xg, yu_g], axis=1)[::-1]   # TE -> LE
    lower = np.stack([xg, yl_g], axis=1)[1:]     # LE+1 -> TE
    return np.vstack([upper, lower])


def build_wing(code='0012', span=6.0, chord=1.0, taper=1.0, sweep_deg=0.0,
               dihedral_deg=0.0, twist_deg=0.0, n_chord=40, n_span=12,
               full=True, te_gap=0.002, cosine=True, section=None):
    """Structured full-wing mesh. x chordwise, y spanwise, z up.

    te_gap: finite trailing-edge thickness (fraction of chord, default
    0.002). Sharp TEs carry a corner singularity that low-order panels
    over-resolve into spurious lift/drag; real airfoils have finite TE
    (~0.2%). Gaps above ~0.003 destabilize the wake coupling; use 0
    only for sharp-TE studies.
    """
    if te_gap > 0.003:
        import warnings
        warnings.warn('te_gap > 0.003 may destabilize the wake coupling')
    if section is not None:
        # custom section loop (TE -> upper -> LE -> lower), e.g. built
        # with section_loop(); replaces the NACA-4 digit section.
        loop = np.asarray(section, float).copy()
    else:
        loop = naca4_section(code, n_chord, cosine=cosine).copy()
    finite_base = te_gap > 0
    if finite_base:
        # target a TOTAL trailing-edge gap of te_gap. Sections that are
        # already open (finite-TE airfoils like DFVLR R4) get closed down;
        # closed sections get opened up. Gaps above ~0.003 destabilize
        # the wake coupling.
        cur = loop[0, 1] - loop[-1, 1]
        adj = (te_gap - cur) / 2.0
        loop[0, 1] += adj
        loop[-1, 1] -= adj
    nl = len(loop)
    sweep = np.tan(np.radians(sweep_deg))
    dih = np.tan(np.radians(dihedral_deg))
    if full:
        etas = np.linspace(-1, 1, 2 * n_span + 1)
    else:
        etas = np.linspace(0, 1, n_span + 1)
    verts = []
    index = {}
    for j, eta in enumerate(etas):
        ay = abs(eta)
        c = chord * (1 - (1 - taper) * ay)
        xoff = sweep * ay * span / 2
        zoff = dih * ay * span / 2
        tw = np.radians(twist_deg * ay)
        R = np.array([[np.cos(tw), -np.sin(tw)], [np.sin(tw), np.cos(tw)]])
        for i, (px, py) in enumerate(loop):
            q = R @ np.array([px * c, py * c])
            verts.append([q[0] + xoff, eta * span / 2, q[1] + zoff])
            index[(i, j)] = len(verts) - 1
    faces = []
    ns = len(etas)
    # NOTE: for sharp TEs the loop closes on itself (loop[-1] == loop[0]),
    # so panels run i -> i+1 without wraparound; the wraparound quad would
    # be a degenerate zero-area panel. A finite te_gap leaves a small open
    # base slot (no base panels: they destabilize conditioning); the leak
    # is negligible for te_gap <= 0.002 and regularizes the TE corner
    # singularity that low-order panels otherwise over-resolve.
    for j in range(ns - 1):
        for i in range(nl - 1):
            faces.append([index[(i, j)], index[(i + 1, j)],
                          index[(i + 1, j + 1)], index[(i, j + 1)]])
    # tip caps (fan around tip loop centroid; no wraparound: TE shared)
    for j in (0, ns - 1):
        pts = np.array([verts[index[(i, j)]] for i in range(nl)])
        c = pts.mean(axis=0)
        ci = len(verts)
        verts.append(c)
        for i in range(nl - 1):
            faces.append([index[(i, j)], index[(i + 1, j)], ci])
    # root cap only for half wing
    mesh = Mesh(np.array(verts), faces)
    # structured topology metadata for wake/strip construction.
    # Wake starts at the UPPER TE corner (index 0): the physical shedding
    # surface at positive lift. (Mid-gap starts detach the wake and break
    # the Kutta coupling; verified by wake-position sensitivity test.)
    te_pts = [np.array(verts[index[(0, j)]]) for j in range(ns)]
    mesh.meta = {'n_loop': nl, 'n_span': ns - 1, 'te_points': te_pts,
                 'span': span, 'chord': chord, 'taper': taper,
                 'sweep_deg': sweep_deg, 'te_gap': te_gap}
    return mesh


def nonplanar_planform(span=1.6, n_st=25, chord_fn=None, xle_fn=None,
                       z_fn=None, twist_fn=None, slope_t_fn=None,
                       slope_c_fn=None, note=''):
    """Span-station description of a general (possibly non-planar) wing.

    This is the geometry INPUT of the TSD Cartesian model: the solver
    rasterises it onto its own (x, y) node grid (TSDSolver.set_planform),
    so taper / sweep / twist / planform shape / per-span section identity
    all enter the field through the same arrays

        xle(y), chord(y), z(y), twist(y)

    and section identity enters as two callables
        slope_t(x_rel, y) -> d(half-thickness)/dx   (local chord units)
        slope_c(x_rel, y) -> d(camber)/dx
    so a thick centrebody and a thin outer panel are one planform.

    Returns dict(y, xle, chord, z, twist_deg, slope_t, slope_c, span, note).
    """
    y = np.linspace(-0.5 * span, 0.5 * span, n_st)
    eta = np.abs(y) / (0.5 * span)
    xle = (np.zeros_like(y) if xle_fn is None
           else np.asarray(xle_fn(eta), float))
    ch = (np.ones_like(y) if chord_fn is None
          else np.asarray(chord_fn(eta), float))
    zz = np.zeros_like(y) if z_fn is None else np.asarray(z_fn(eta), float)
    tw = (np.zeros_like(y) if twist_fn is None
          else np.asarray(twist_fn(eta), float))
    return {'y': y, 'xle': xle, 'chord': ch, 'z': zz, 'twist_deg': tw,
            'slope_t': slope_t_fn, 'slope_c': slope_c_fn, 'span': span,
            'eta': eta, 'note': note}


def thickness_slope_np(x, t=0.12):
    """NACA-4 half-thickness slope, vectorised over a thickness array."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    t = np.asarray(t, dtype=float)
    return 5 * t * (0.2969 / (2 * np.sqrt(np.maximum(x, 1e-9))) - 0.1260
                    - 2 * 0.3516 * x + 3 * 0.2843 * x ** 2
                    - 4 * 0.1036 * x ** 3)


def flying_wing_planform(span=3.0, chord_root=1.0, taper=0.18, eta_crank=0.22,
                         sweep_le_deg=32.0, sweep_le_in_deg=8.0,
                         dihedral_deg=0.0, washout_deg=3.0, n_st=41,
                         t_centre=0.16, t_outer=0.10):
    """A REPRESENTATIVE flying-wing planform (B-2-like), NOT a B-2.

    No accurate public B-2 geometry exists, so this is a parametric
    blended flying wing built to exercise every capability the planar
    model lacked: a swept leading edge with an inboard crank, a blended
    centrebody (constant chord inboard), strong taper to a small tip
    chord, washout, optional dihedral, and a THICKER centre section than
    the outer panels (section identity varying with span).

    Numbers chosen: AR ~ 3 (a low-AR blended body), LE sweep 32 deg
    outboard / 8 deg through the centrebody, tip chord 18% of root,
    centre section 16% thick vs 10% outboard. ALL REPRESENTATIVE.
    """
    def chord_fn(eta):
        e = np.clip((eta - eta_crank) / max(1.0 - eta_crank, 1e-9), 0, 1)
        return chord_root * (1.0 - (1.0 - taper) * e)

    def xle_fn(eta):
        yv = eta * 0.5 * span
        yc = eta_crank * 0.5 * span
        return np.where(
            yv <= yc, np.tan(np.radians(sweep_le_in_deg)) * yv,
            np.tan(np.radians(sweep_le_in_deg)) * yc
            + np.tan(np.radians(sweep_le_deg)) * (yv - yc))

    def z_fn(eta):
        return np.tan(np.radians(dihedral_deg)) * eta * 0.5 * span

    def twist_fn(eta):
        return -washout_deg * eta

    def slope_t_fn(xrel, yy):
        eta = np.abs(np.asarray(yy, float)) / (0.5 * span)
        e = np.clip((eta - eta_crank) / max(1.0 - eta_crank, 1e-9), 0, 1)
        t = t_centre + (t_outer - t_centre) * e
        return thickness_slope_np(xrel, t)

    def slope_c_fn(xrel, yy):
        return np.zeros_like(np.asarray(xrel, float))

    return nonplanar_planform(
        span=span, n_st=n_st, chord_fn=chord_fn, xle_fn=xle_fn,
        z_fn=z_fn, twist_fn=twist_fn, slope_t_fn=slope_t_fn,
        slope_c_fn=slope_c_fn,
        note='representative flying wing (not a B-2)')


def read_stl(path):
    """Read ASCII or binary STL -> Mesh (triangles, welded)."""
    with open(path, 'rb') as fh:
        data = fh.read()
    tris = None
    if len(data) > 84:
        count = int.from_bytes(data[80:84], 'little')
        if len(data) == 84 + count * 50:
            raw = np.frombuffer(data[84:], dtype=np.uint8).reshape(count, 50)
            f32 = raw[:, :48].view(np.float32)
            tris = f32[:, 3:12].reshape(count, 3, 3)
    if tris is None:
        text = data.decode('utf-8', errors='ignore')
        import re
        nums = re.findall(r'vertex\s+([^\n]+)', text)
        pts = np.array([[float(x) for x in s.split()] for s in nums])
        tris = pts.reshape(-1, 3, 3)
    # weld
    key = np.round(tris.reshape(-1, 3), 9)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    verts = np.unique(key, axis=0)
    faces = inv.reshape(-1, 3).tolist()
    # drop degenerate
    good = []
    for f in faces:
        if len(set(f)) == 3:
            good.append(f)
    return Mesh(verts, good)


def decimate(mesh, target=2000, seed=0):
    """Clustering decimation to ~target panels. Returns polygon Mesh."""
    rng = np.random.default_rng(seed)
    c = mesh.centroid
    lo, hi = c.min(axis=0), c.max(axis=0)
    ext = float(np.linalg.norm(mesh.area.sum()))
    ncell = max(1, int(round((mesh.npanels / target) ** (1 / 3))))
    # k-means-lite: uniform grid buckets, split biggest until target
    buckets = [np.arange(mesh.npanels)]
    while len(buckets) < target:
        sizes = np.array([len(b) for b in buckets])
        i = int(np.argmax(sizes))
        if sizes[i] < 2:
            break
        b = buckets.pop(i)
        spread = c[b].max(axis=0) - c[b].min(axis=0)
        ax = int(np.argmax(spread))
        order = np.argsort(c[b][:, ax])
        mid = len(order) // 2
        buckets.append(b[order[:mid]])
        buckets.append(b[order[mid:]])
    verts = []
    faces = []
    for b in buckets:
        # refit plane by SVD on member vertices
        idx = np.unique([v for f in [mesh.faces[k] for k in b] for v in f])
        p = mesh.verts[idx]
        ctr = p.mean(axis=0)
        _, _, vt = np.linalg.svd(p - ctr, full_matrices=False)
        n = vt[-1]
        if n @ mesh.normal[b].sum(axis=0) < 0:
            n = -n
        u, w = vt[0], vt[1]
        q = (p - ctr)
        ang = np.arctan2(q @ w, q @ u)
        hull = np.argsort(ang)
        ring = []
        for k in hull:
            ring.append(len(verts))
            verts.append(ctr + (q[k] @ u) * u + (q[k] @ w) * w)
        faces.append(ring)
    return Mesh(np.array(verts), faces)
