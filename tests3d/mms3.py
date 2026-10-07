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

S0 = np.array([2.0, 0.0, 0.0])
Q = 5.0


def phi_star(P):
    P = np.asarray(P, dtype=float)
    return P[..., 0] + Q / (4 * np.pi * np.linalg.norm(P - S0, axis=-1))


def grad_star(P):
    P = np.asarray(P, dtype=float)
    r = P - S0
    d = np.linalg.norm(r, axis=-1, keepdims=True)
    return np.array([1., 0., 0.]) - Q * r / (4 * np.pi * d ** 3)


# exact jumps for TRUE representation Phi = Phi_inf + B.sigma + C.mu,
# V = Vinf + (negated kernel sums) [current code convention]:
#   sigma: single-layer jump [dPhi/dn] = -sigma  =>  sigma = Vin.n - dPhi*_out/dn
#   mu: doublet jump [Phi] = +mu  =>  mu = Phi*_out - Phi_inf
sigma = (N @ np.array([1., 0., 0.])) - np.einsum('ij,ij->i', grad_star(C), N)
mu = phi_star(C) - C[:, 0]

pc = _precompute(m)
size = np.sqrt(m.area)
print('== exterior potential ==')
for P in [np.array([3., 0, 0]), np.array([0., 0, 2.]),
          np.array([-1.5, 0.5, 0.3]), np.array([1.3, 0.0, 0.0])]:
    om = _solid_angle_rows(P, pc['fan'], pc['fmask'])
    b = _source_pot_rows(P, pc['qp'], pc['qw'], m.area, size)
    num = P[0] + b.dot(sigma) + ((-om / (4 * np.pi)).dot(mu))
    print('P', P, 'num %.4f exact %.4f' % (num, phi_star(P)))

print('== surface velocity ==')
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
g = grad_star(p)
gt = g - (g @ N[i]) * N[i]
print('kernel |Vt| =', round(float(np.linalg.norm(v)), 4),
      ' exact |Vt| =', round(float(np.linalg.norm(gt)), 4))
