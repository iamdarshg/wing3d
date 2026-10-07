import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh
from wing3d.solver import (assemble, build_wake, _precompute,
                            _source_vel_rows, _loop_vel_rows,
                            wake_precompute)

nc, ns = 20, 12
xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, nc + 1)))
ys = np.linspace(-3, 3, ns + 1)
verts = []
up_id, lo_id = {}, {}
for j, y in enumerate(ys):
    for i, x in enumerate(xs):
        up_id[(i, j)] = len(verts)
        verts.append([x, y, 0.0])
for j, y in enumerate(ys):
    for i, x in enumerate(xs):
        lo_id[(i, j)] = len(verts)
        verts.append([x, y, -1e-6])
faces = []
for j in range(ns):
    for i in range(nc):
        faces.append([up_id[(i, j)], up_id[(i + 1, j)],
                      up_id[(i + 1, j + 1)], up_id[(i, j + 1)]])
for j in range(ns):
    for i in range(nc):
        faces.append([lo_id[(i, j)], lo_id[(i, j + 1)],
                      lo_id[(i + 1, j + 1)], lo_id[(i + 1, j)]])
m = Mesh(np.array(verts), faces)
nup = ns * nc
info = []
for j in range(ns):
    up = j * nc + (nc - 1)
    lo = nup + j * nc + (nc - 1)
    info.append((up, lo, np.array([1.0, ys[j], 0.0]),
                 np.array([1.0, ys[j + 1], 0.0])))
wakes = build_wake(m, info)
Afull, Brow = assemble(m, wakes, kutta_mode='doublet')
n = m.npanels
s = len(wakes)
Abody = Afull[:n, :n]
Wcol = Afull[:n, n:]
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
sigma = (m.normal @ V)
rhs0 = (Brow @ sigma)  # current code: RHS=+Bsigma
pc = _precompute(m)
wc = wake_precompute(wakes)
size = np.sqrt(m.area)


def CL_of(mu, muw):
    vel = np.zeros((n, 3))
    for i in range(n):
        p = m.centroid[i]
        v = V - (sigma[:, None] * _source_vel_rows(
            p, pc['qp'], pc['qw'], m.area, size,
            self_idx=i, normal=m.normal[i])).sum(axis=0) \
            - (mu[:, None] * _loop_vel_rows(
                p, pc['ed'], pc['emask'])).sum(axis=0) \
            - (muw[wc['owner'], None] * _loop_vel_rows(
                p, wc['ed'], wc['emask'])).sum(axis=0)
        v = v - 0.5 * sigma[i] * m.normal[i]
        vel[i] = v
    vel = vel - (np.einsum('ij,ij->i', vel, m.normal))[:, None] * m.normal
    cp = 1 - np.einsum('ij,ij->i', vel, vel)
    return float((-cp[:, None] * m.normal * m.area[:, None]).sum(axis=0)[2] / 6.0)


mu0 = np.linalg.solve(Abody, rhs0)
for dmw in [0.0, 0.2, -0.2]:
    muw = np.full(s, dmw)
    mu = np.linalg.solve(Abody, rhs0 - Wcol @ muw)
    print('prescribed muw=%+.2f CL=%+.4f' % (dmw, CL_of(mu, muw)))
