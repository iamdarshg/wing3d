import numpy as np
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
mu = np.linalg.solve(A, -(Brow @ sigma))
pc = _precompute(m)
size = np.sqrt(m.area)
i = int(np.argmin(abs(m.centroid[:, 0])))
p = m.centroid[i]
n = m.normal[i]
t = V - (V @ n) * n
t /= np.linalg.norm(t)


def brute_vel(j, N=30):
    """True velocity (physical, outward-positive source) from panel j."""
    v = np.asarray(m.verts[m.faces[j]])
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
                    rv = p - g
                    dd = np.linalg.norm(rv)
                    if dd < 1e-9:
                        continue
                    # physical source velocity: +sigma*A*rhat/4pir^2
                    out += ar * rv / dd ** 3
                if ii + jj + 1 < N:
                    Q0 = a0 + (b0 - a0) * (ii + 1) / N + (c0 - a0) * jj / N
                    Q1 = a0 + (b0 - a0) * ii / N + (c0 - a0) * (jj + 1) / N
                    Q2 = a0 + (b0 - a0) * (ii + 1) / N + (c0 - a0) * (jj + 1) / N
                    g = (Q0 + Q1 + Q2) / 3
                    ar = 0.5 * np.linalg.norm(np.cross(Q1 - Q0, Q2 - Q0))
                    rv = p - g
                    dd = np.linalg.norm(rv)
                    if dd < 1e-9:
                        continue
                    out += ar * rv / dd ** 3
    return out / (4 * np.pi)


K = _source_vel_rows(p, pc['qp'], pc['qw'], m.area, size)
d = np.linalg.norm(m.centroid - p, axis=1)
print('tangential +physical-source-velocity per panel: brute vs kernel(-sign)')
for j in np.argsort(d)[1:7]:
    b = float(brute_vel(j) @ t) * sigma[j]
    k = float(K[j] @ t) * sigma[j]
    print(f'panel {j} dist={d[j]:.3f} brute={b:+.5f} kernel={k:+.5f} ratio={b / k if abs(k) > 1e-9 else float("nan"):.2f}')
