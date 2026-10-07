import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import (_precompute, _solid_angle_rows, _source_pot_rows,
                            _source_vel_rows, _loop_vel_rows)

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
C = m.centroid
N = m.normal
ct = C[:, 0]  # cos theta (unit sphere)
sigma = ct
mu = 0.5 * ct
pc = _precompute(m)
size = np.sqrt(m.area)

print('== exterior potential ==')
for P in [np.array([2., 0, 0]), np.array([0., 0, 2.]),
          np.array([1.2, 0.3, -0.2])]:
    om = _solid_angle_rows(P, pc['fan'], pc['fmask'])
    b = _source_pot_rows(P, pc['qp'], pc['qw'], m.area, size)
    num = P[0] + b.dot(sigma) + ((-om / (4 * np.pi)).dot(mu))
    r = np.linalg.norm(P)
    exact = (r + 0.5 / r ** 2) * (P[0] / r)
    print('P', P, 'num %.4f exact %.4f' % (num, exact))

print('== surface velocity (negated assembly) ==')
i = int(np.argmin(abs(C[:, 0])))
p = C[i]
v = (np.array([1., 0., 0.])
     - (sigma[:, None] * _source_vel_rows(
         p, pc['qp'], pc['qw'], m.area, size,
         self_idx=i, normal=N[i])).sum(axis=0)
     - (mu[:, None] * _loop_vel_rows(p, pc['ed'], pc['emask'])).sum(axis=0))
v = v - (v @ N[i]) * N[i]
t = np.array([1., 0., 0.])
t = t - (t @ N[i]) * N[i]
t /= np.linalg.norm(t)
print('|Vt| =', round(float(np.linalg.norm(v)), 3), '(analytic 1.5*sin67=1.38)')
