import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import (_precompute, _source_vel_rows, assemble)

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
i = int(np.argmin(abs(m.centroid[:, 0])))
p = m.centroid[i]
n = m.normal[i]
t = V - (V @ n) * n
t /= np.linalg.norm(t)


def brute_panel(j, N=36):
    v = np.asarray(m.faces[j] and m.verts[m.faces[j]])
    a0 = v[0]
    out = np.zeros(3)
    for k in range(1, len(v) - 1):
        b0, c0 = v[k], v[k + 1]
        for ii in range(N):
            for jj in range(N - ii):
                for (P0, P1, P2) in [((ii, jj), (ii + 1, jj), (ii, jj + 1))]:
                    q0 = a0 + (b0 - a0) * P0[0] / N + (c0 - a0) * P0[1] / N
                    q1 = a0 + (b0 - a0) * P1[0] / N + (c0 - a0) * P1[1] / N
                    q2 = a0 + (b0 - a0) * P2[0] / N + (c0 - a0) * P2[1] / N
                    g = (q0 + q1 + q2) / 3
                    ar = 0.5 * np.linalg.norm(np.cross(q1 - q0, q2 - q0))
                    rv = g - p
                    dd = np.linalg.norm(rv)
                    if dd < 1e-9:
                        continue
                    out += ar * rv / dd ** 3
                if ii + jj + 1 < N:
                    Q0 = a0 + (b0 - a0) * (ii + 1) / N + (c0 - a0) * jj / N
                    Q1 = a0 + (b0 - a0) * ii / N + (c0 - a0) * (jj + 1) / N
                    Q2 = a0 + (b0 - a0) * (ii + 1) / N + (c0 - a0) * (jj + 1) / N
                    g = (Q0 + Q1 + Q2) / 3
                    ar = 0.5 * np.linalg.norm(np.cross(Q1 - Q0, Q2 - Q0))
                    rv = g - p
                    dd = np.linalg.norm(rv)
                    if dd < 1e-9:
                        continue
                    out += ar * rv / dd ** 3
    return out / (4 * np.pi)


K = _source_vel_rows(p, pc['qp'], pc['qw'], m.area, size)
d = np.linalg.norm(m.centroid - p, axis=1)
near = (d < 0.6) & (np.arange(len(d)) != i)
b_sum = np.zeros(3)
for j in np.where(near)[0]:
    b_sum += brute_panel(j) * sigma[j]
k_sum = (sigma[:, None] * K)[near].sum(axis=0)
print('near-field source .t: brute %.5f kernel %.5f' % (b_sum @ t, k_sum @ t))
print('n near panels:', near.sum())
