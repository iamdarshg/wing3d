import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import (_precompute, _solid_angle_rows, _source_pot_rows)

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

S0 = np.array([5.0, 0.0, 0.0])
Q = 2.0
C = m.centroid
N = m.normal


def phi_star(P):
    P = np.asarray(P, dtype=float)
    return P[..., 0] + Q / (4 * np.pi * np.linalg.norm(P - S0, axis=-1))


def grad_star(P):
    P = np.asarray(P, dtype=float)
    r = P - S0
    d = np.linalg.norm(r, axis=-1, keepdims=True)
    return np.array([1., 0., 0.]) - Q * r / (4 * np.pi * d ** 3)


# exact jumps for uniform-interior representation (my C convention:
# doublet jump [Phi] = +mu, single-layer [dPhi/dn] = -sigma)
sigma = (N @ np.array([1., 0., 0.])) - np.einsum('ij,ij->i', grad_star(C), N)
mu = phi_star(C) - C[:, 0]
print('sigma range:', sigma.min().round(3), sigma.max().round(3))
print('mu range:', mu.min().round(3), mu.max().round(3))

pc = _precompute(m)
size = np.sqrt(m.area)
for P in [np.array([2., 0, 0]), np.array([0., 0, 2.]),
          np.array([-1.5, 0.5, 0.3]), np.array([0., 1.5, -0.5]),
          np.array([1.1, 0.2, 0.1])]:
    om = _solid_angle_rows(P, pc['fan'], pc['fmask'])
    b = _source_pot_rows(P, pc['qp'], pc['qw'], m.area, size)
    num = P[0] + b.dot(sigma) + ((-om / (4 * np.pi)).dot(mu))
    print('P', P, 'num %.4f exact %.4f' % (num, phi_star(P)))
