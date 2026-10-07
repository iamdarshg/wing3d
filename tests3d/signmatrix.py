import numpy as np
import sys
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
sigma = (m.normal @ V)
x = np.linalg.solve(A, -(Brow @ sigma))
mu = x
pc = _precompute(m)
size = np.sqrt(m.area)
n = m.npanels
C = m.centroid


def velocities(sgn, self_mode):
    vel = np.zeros((n, 3))
    for i in range(n):
        p = C[i]
        sv = (sigma[:, None] * _source_vel_rows(
            p, pc['qp'], pc['qw'], m.area, size,
            self_idx=i, normal=m.normal[i])).sum(axis=0)
        lv = (mu[:, None] * _loop_vel_rows(p, pc['ed'], pc['emask'])).sum(axis=0)
        v = V + sgn * (sv + lv)
        if self_mode == 'plus':
            v += 0.5 * sigma[i] * m.normal[i]
        elif self_mode == 'minus':
            v -= 0.5 * sigma[i] * m.normal[i]
        vel[i] = v
    return vel


th = np.arccos(np.clip(C[:, 0], -1, 1))
exact = 1 - 2.25 * np.sin(th) ** 2
# NOTE _source_vel_rows currently zeroes the self row; modes add it back
for sgn in [+1, -1]:
    for mode in ['plus', 'minus', 'zero']:
        vel = velocities(sgn, mode)
        vel = vel - (np.einsum('ij,ij->i', vel, m.normal))[:, None] * m.normal
        cp = 1 - np.einsum('ij,ij->i', vel, vel)
        rms = float(np.sqrt(((cp - exact) ** 2).mean()))
        print(f'sgn={sgn:+d} self={mode:5s} rms={rms:.4f}')
