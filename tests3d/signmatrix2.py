import numpy as np
from wing3d.geometry import Mesh
from wing3d.solver import (_precompute, _loop_vel_rows, assemble,
                            _centroids)

nlat, nlon = 12, 24
verts = []
for i in range(nlat + 1):
    th = np.pi * i / nlat
    for j in range(nlon):
        ph = 2 * np.pi * j / nlon
        verts.append([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                      np.cos(th)])
verts = np.array(verts)
faces = []
for i in range(nlat):
    for j in range(nlon):
        a = i * nlon + j
        b = i * nlon + (j + 1) % nlon
        c = (i + 1) * nlon + j
        d = (i + 1) * nlon + (j + 1) % nlon
        if i > 0:
            faces.append([a, c, b])
        if i < nlat - 1:
            faces.append([b, c, d])
m = Mesh(verts, faces)
A, Brow = assemble(m, [])
V = np.array([1., 0., 0.])
sigma = (m.normal @ V)
mu = np.linalg.solve(A, -(Brow @ sigma))
pc = _precompute(m)
size = np.sqrt(m.area)
n = m.npanels
C = m.centroid
th = np.arccos(np.clip(C[:, 0], -1, 1))
exact = 1 - 2.25 * np.sin(th) ** 2


def src_vel(Pp, i, far_sign, near_sign):
    c = _centroids(pc['qp'], pc['qw'])
    r = np.linalg.norm(c - Pp, axis=1)
    far = r > 4 * size
    far[i] = False
    out = np.zeros((len(m.area), 3))
    rv = Pp - c[far]
    d = np.maximum(r[far], 1e-14)
    out[far] = far_sign * m.area[far][:, None] * rv / (4 * np.pi * d[:, None] ** 3)
    near = ~far
    near[i] = False
    dd = pc['qp'][near] - Pp
    dist = np.maximum(np.linalg.norm(dd, axis=-1), 1e-12)
    out[near] = near_sign * (-np.sum((pc['qw'][near] / dist ** 3)[..., None] * dd,
                                     axis=1) / (4 * np.pi))
    return out


for fs in [+1, -1]:
    for ns_ in [+1, -1]:
        for ds in [+1, -1]:
            vel = np.zeros((n, 3))
            for i in range(n):
                p = C[i]
                sv = (sigma[:, None] * src_vel(p, i, fs, ns_)).sum(axis=0)
                lv = (mu[:, None] * _loop_vel_rows(
                    p, pc['ed'], pc['emask'])).sum(axis=0)
                vel[i] = V + sv + ds * lv
            vel = vel - (np.einsum('ij,ij->i', vel, m.normal))[:, None] * m.normal
            cp = 1 - np.einsum('ij,ij->i', vel, vel)
            print(f'far={fs:+d} near={ns_:+d} loop={ds:+d} '
                  f'rms={np.sqrt(((cp - exact) ** 2).mean()):.4f}')
