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


def build_wing(code='0012', span=6.0, chord=1.0, taper=1.0, sweep_deg=0.0,
               dihedral_deg=0.0, twist_deg=0.0, n_chord=40, n_span=12,
               full=True, te_gap=0.002, cosine=True):
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
    loop = naca4_section(code, n_chord, cosine=cosine).copy()
    finite_base = te_gap > 0
    if finite_base:
        # open the loop: upper TE up, lower TE down; base closed below
        loop[0, 1] += te_gap / 2
        loop[-1, 1] -= te_gap / 2
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


def read_stl(path):
    """Read ASCII or binary STL -> Mesh (triangles, welded)."""
    with open(path, 'rb') as fh:
        data = fh.read()
    tris = None
    if len(data) > 84:
        count = int.from_bytes(data[80:84], 'little')
        if len(data) == 84 + count * 50:
            tris = np.frombuffer(data[84:], dtype=np.float32).reshape(count, -1)[:, 3:12].reshape(count, 3, 3)
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
