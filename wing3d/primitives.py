"""Primitive 3D meshes (quads/tris) + merge/weld. All output closed Mesh objects."""
import numpy as np
from .geometry import Mesh


def _grid(nu, nv):
    return np.meshgrid(np.linspace(0, 1, nu + 1), np.linspace(0, 1, nv + 1),
                       indexing='ij')


def box(size=(1.0, 1.0, 1.0), n=(8, 8, 8), center=(0.0, 0.0, 0.0)):
    """Subdivided box. size=(sx,sy,sz), n=(nx,ny,nz) panels per face pair."""
    sx, sy, sz = size
    nx, ny, nz = n
    cx, cy, cz = center
    verts, faces = [], []

    def add_patch(corner, eu, ev, nu, nv):
        base = len(verts)
        U, V = _grid(nu, nv)
        for i in range(nu + 1):
            for j in range(nv + 1):
                verts.append(corner + eu * U[i, j] + ev * V[i, j])
        for i in range(nu):
            for j in range(nv):
                a = base + i * (nv + 1) + j
                faces.append([a, a + nv + 1, a + nv + 2, a + 1])

    x0, y0, z0 = cx - sx / 2, cy - sy / 2, cz - sz / 2
    # +x, -x, +y, -y, +z, -z (outward CCW ensured by Mesh auto-orient)
    add_patch(np.array([x0 + sx, y0, z0]), np.array([0, sy, 0]),
              np.array([0, 0, sz]), ny, nz)
    add_patch(np.array([x0, y0, z0]), np.array([0, 0, sz]),
              np.array([0, sy, 0]), nz, ny)
    add_patch(np.array([x0, y0 + sy, z0]), np.array([sx, 0, 0]),
              np.array([0, 0, sz]), nx, nz)
    add_patch(np.array([x0, y0, z0]), np.array([0, 0, sz]),
              np.array([sx, 0, 0]), nz, nx)
    add_patch(np.array([x0, y0, z0 + sz]), np.array([sx, 0, 0]),
              np.array([0, sy, 0]), nx, ny)
    add_patch(np.array([x0, y0, z0]), np.array([0, sy, 0]),
              np.array([sx, 0, 0]), ny, nx)
    return Mesh(np.array(verts), faces)


def ellipsoid(a=1.0, b=0.5, c=0.5, nlat=16, nlon=32, center=(0, 0, 0)):
    """Axis-aligned ellipsoid (clean quads, pole fans)."""
    verts = []
    for i in range(nlat + 1):
        th = np.pi * i / nlat
        for j in range(nlon):
            ph = 2 * np.pi * j / nlon
            verts.append([a * np.sin(th) * np.cos(ph) + center[0],
                          b * np.sin(th) * np.sin(ph) + center[1],
                          c * np.cos(th) + center[2]])
    verts = np.array(verts)
    faces = []
    for i in range(nlat):
        for j in range(nlon):
            a_ = i * nlon + j
            b_ = i * nlon + (j + 1) % nlon
            c_ = (i + 1) * nlon + j
            d_ = (i + 1) * nlon + (j + 1) % nlon
            if i > 0:
                faces.append([a_, c_, b_])
            if i < nlat - 1:
                faces.append([b_, c_, d_])
    # pole caps are the i=0/i=nlat-1 triangle rings above (degenerate quads
    # avoided by using triangles at poles)
    return Mesh(verts, faces)


def cylinder(radius=0.5, length=2.0, nrad=24, nlen=12, center=(0, 0, 0),
             axis='x'):
    """Capped cylinder along axis."""
    verts = []
    for j in range(nlen + 1):
        s = -length / 2 + length * j / nlen
        for i in range(nrad):
            ph = 2 * np.pi * i / nrad
            p = [radius * np.cos(ph), radius * np.sin(ph), s]
            if axis == 'x':
                p = [s, radius * np.cos(ph), radius * np.sin(ph)]
            elif axis == 'y':
                p = [radius * np.cos(ph), s, radius * np.sin(ph)]
            verts.append(np.array(p) + np.array(center))
    verts = np.array(verts)
    faces = []
    for j in range(nlen):
        for i in range(nrad):
            i2 = (i + 1) % nrad
            faces.append([j * nrad + i, j * nrad + i2,
                          (j + 1) * nrad + i2, (j + 1) * nrad + i])
    for j, sgn in ((0, -1), (nlen, 1)):
        ring = [j * nrad + i for i in range(nrad)]
        ctr = np.mean(verts[ring], axis=0)
        ci = len(verts)
        verts = np.vstack([verts, ctr])
        for i in range(nrad):
            i2 = (i + 1) % nrad
            faces.append([ring[i], ring[i2], ci])
    return Mesh(verts, faces)


def frustum(r0=0.5, r1=0.2, length=1.0, nrad=24, nlen=6, center=(0, 0, 0),
            axis='x', cap_back=True):
    """Cone frustum (nose cone with r1 front). Open front, optional back cap."""
    verts = []
    for j in range(nlen + 1):
        s = -length / 2 + length * j / nlen
        r = r0 + (r1 - r0) * j / nlen
        for i in range(nrad):
            ph = 2 * np.pi * i / nrad
            p = [r * np.cos(ph), r * np.sin(ph), s]
            if axis == 'x':
                p = [s, r * np.cos(ph), r * np.sin(ph)]
            elif axis == 'y':
                p = [r * np.cos(ph), s, r * np.sin(ph)]
            verts.append(np.array(p) + np.array(center))
    verts = np.array(verts)
    faces = []
    for j in range(nlen):
        for i in range(nrad):
            i2 = (i + 1) % nrad
            faces.append([j * nrad + i, j * nrad + i2,
                          (j + 1) * nrad + i2, (j + 1) * nrad + i])
    if cap_back:
        ring = [(nlen) * nrad + i for i in range(nrad)]
        ctr = np.mean(verts[ring], axis=0)
        ci = len(verts)
        verts = np.vstack([verts, ctr])
        for i in range(nrad):
            i2 = (i + 1) % nrad
            faces.append([ring[i], ring[i2], ci])
    return Mesh(verts, faces)


def wedge(x0=0.0, x1=1.0, y0=-0.5, y1=0.5, z0=0.0, z1=0.3, nx=8, ny=8):
    """Ramp wedge (Ahmed slantback building block): footprint [x0,x1]x[y0,y1],
    height z0 at x0 rising to z1 at x1. Closed with bottom + sides."""
    verts, faces = [], []

    def sheet(xa, xb, za, zb, nu):
        base = len(verts)
        for i in range(nu + 1):
            x = xa + (xb - xa) * i / nu
            z = za + (zb - za) * i / nu
            for j in range(ny + 1):
                y = y0 + (y1 - y0) * j / ny
                verts.append([x, y, z])
        for i in range(nu):
            for j in range(ny):
                a = base + i * (ny + 1) + j
                faces.append([a, a + ny + 1, a + ny + 2, a + 1])
        return base

    sheet(x0, x1, z0, z1, nx)  # slanted top
    # bottom
    b = len(verts)
    for i in range(nx + 1):
        x = x0 + (x1 - x0) * i / nx
        for j in range(ny + 1):
            y = y0 + (y1 - y0) * j / ny
            verts.append([x, y, z0 - 0.001])
    for i in range(nx):
        for j in range(ny):
            a = b + i * (ny + 1) + j
            faces.append([a, a + 1, a + ny + 2, a + ny + 1])
    # four side walls (winding: outward normals; verified by watertight
    # + outward-normal check -- inward faces break Morino globally)
    # front wall x=x0 (outward -x)
    b = len(verts)
    for j in range(ny + 1):
        y = y0 + (y1 - y0) * j / ny
        verts.append([x0, y, z0 - 0.001])
        verts.append([x0, y, z0])
    for j in range(ny):
        a = b + 2 * j
        faces.append([a, a + 1, a + 3, a + 2])
    # back wall x=x1 (outward +x)
    b = len(verts)
    for j in range(ny + 1):
        y = y0 + (y1 - y0) * j / ny
        verts.append([x1, y, z0 - 0.001])
        verts.append([x1, y, z1])
    for j in range(ny):
        a = b + 2 * j
        faces.append([a, a + 2, a + 3, a + 1])
    # side walls y=y0 (outward -y), y=y1 (outward +y)
    b = len(verts)
    for i in range(nx + 1):
        x = x0 + (x1 - x0) * i / nx
        z = z0 + (z1 - z0) * i / nx
        verts.append([x, y0, z0 - 0.001])
        verts.append([x, y0, z])
    for i in range(nx):
        a = b + 2 * i
        faces.append([a, a + 2, a + 3, a + 1])
    b = len(verts)
    for i in range(nx + 1):
        x = x0 + (x1 - x0) * i / nx
        z = z0 + (z1 - z0) * i / nx
        verts.append([x, y1, z0 - 0.001])
        verts.append([x, y1, z])
    for i in range(nx):
        a = b + 2 * i
        faces.append([a, a + 1, a + 3, a + 2])
    return Mesh(np.array(verts), faces)


def transform(mesh, scale=1.0, translate=(0, 0, 0), rot_z_deg=0.0):
    """Scale + translate (+ optional z-rotation). Returns new Mesh."""
    v = np.asarray(mesh.verts, dtype=float) * scale
    if rot_z_deg:
        t = np.radians(rot_z_deg)
        R = np.array([[np.cos(t), -np.sin(t), 0],
                      [np.sin(t), np.cos(t), 0], [0, 0, 1]])
        v = v @ R.T
    v = v + np.asarray(translate, dtype=float)
    return Mesh(v, [list(f) for f in mesh.faces])


def merge(meshes):
    """Concatenate meshes into one (welds coincident verts for neighbor search)."""
    verts, faces = [], []
    for m in meshes:
        base = len(verts)
        verts.extend([tuple(p) for p in np.asarray(m.verts)])
        for f in m.faces:
            faces.append([base + k for k in f])
    verts = np.array(verts)
    # weld duplicates (rounded) so neighbor search sees shared verts
    key = np.round(verts, 9)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    faces = [[int(inv[k]) for k in f] for f in faces]
    # drop degenerate
    good = []
    for f in faces:
        if len(set(f)) >= 3:
            good.append(f)
    return Mesh(np.unique(key, axis=0), good)
