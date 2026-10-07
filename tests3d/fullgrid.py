import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import (_precompute, _source_vel_rows, _loop_vel_rows,
                            assemble)

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
pc = _precompute(m)
size = np.sqrt(m.area)
n = m.npanels
C = m.centroid
th = np.arccos(np.clip(C[:, 0], -1, 1))
exact = 1 - 2.25 * np.sin(th) ** 2
for ssign in [+1, -1]:
    sigma = ssign * (m.normal @ V)
    for rsign in [-1, +1]:
        mu = np.linalg.solve(A, rsign * (Brow @ sigma))
        for asign in [+1, -1]:
            vel = np.zeros((n, 3))
            for i in range(n):
                p = C[i]
                v = V + asign * (
                    (sigma[:, None] * _source_vel_rows(
                        p, pc['qp'], pc['qw'], m.area, size,
                        self_idx=i, normal=m.normal[i])).sum(axis=0)
                    + (mu[:, None] * _loop_vel_rows(
                        p, pc['ed'], pc['emask'])).sum(axis=0))
                vel[i] = v
            vel = (vel - (np.einsum('ij,ij->i', vel, m.normal))[:, None]
                   * m.normal)
            cp = 1 - np.einsum('ij,ij->i', vel, vel)
            print(f'sig={ssign:+d} RHS={rsign:+d} ASM={asign:+d} '
                  f'rms={np.sqrt(((cp - exact) ** 2).mean()):.4f}')
