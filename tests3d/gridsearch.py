import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
import wing3d.solver as S
from wing3d.geometry import Mesh
from wing3d.solver import assemble, build_wake

nc, ns = 12, 8
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
    info.append((j * nc + (nc - 1), nup + j * nc + (nc - 1),
                 np.array([1.0, ys[j], 0.0]), np.array([1.0, ys[j + 1], 0.0])))
wakes = build_wake(m, info)
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
n = m.npanels
s = len(wakes)
pc = S._precompute(m)
wc = S.wake_precompute(wakes)
size = np.sqrt(m.area)
C = m.centroid


def run(rhs_sign, asm_sign, ksign, kmode):
    A, Brow = assemble(m, wakes, kutta_mode=kmode, kutta_sign=ksign)
    sigma = (m.normal @ V)
    rhs = np.zeros(n + s)
    rhs[:n] = rhs_sign * (Brow @ sigma)
    if kmode == 'potential':
        kr, _ = S.kutta_rhs(m, wakes, V, Brow)
        # kutta_rhs currently returns new-theory values; adjust by theory below
        rhs[n:] = kr
    x = np.linalg.solve(A, rhs)
    mu, muw = x[:n], x[n:]
    vel = np.zeros((n, 3))
    for i in range(n):
        p = C[i]
        v = V + asm_sign * (
            (sigma[:, None] * S._source_vel_rows(
                p, pc['qp'], pc['qw'], m.area, size,
                self_idx=i, normal=m.normal[i])).sum(axis=0)
            + (mu[:, None] * S._loop_vel_rows(
                p, pc['ed'], pc['emask'])).sum(axis=0)
            + (muw[wc['owner'], None] * S._loop_vel_rows(
                p, wc['ed'], wc['emask'])).sum(axis=0))
        vel[i] = v
    vel = vel - (np.einsum('ij,ij->i', vel, m.normal))[:, None] * m.normal
    cp = 1 - np.einsum('ij,ij->i', vel, vel)
    F = (-cp[:, None] * m.normal * m.area[:, None]).sum(axis=0)
    return float(F[2] / 3.0), float(F[0] / 3.0)


print('(thin-airfoil AR6 alpha=4 wants CL=+0.33, CDp~0)')
for kmode in ['doublet']:
    for rhs_sign in [+1, -1]:
        for asm_sign in [+1, -1]:
            for ksign in [+1, -1]:
                try:
                    cl, cdp = run(rhs_sign, asm_sign, ksign, kmode)
                    print(f'{kmode} RHS={rhs_sign:+d} ASM={asm_sign:+d} K={ksign:+d} CL={cl:+.3f} CDp={cdp:+.4f}')
                except Exception as e:
                    print(kmode, rhs_sign, asm_sign, ksign, 'FAIL', e)
