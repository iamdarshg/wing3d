import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
import wing3d.solver as S
from wing3d.geometry import build_wing
from wing3d.solver import assemble, build_wake_from_meta
from wing3d.forces import pressure_forces

mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=16, n_span=6)
wakes = build_wake_from_meta(mesh)
a = np.radians(4)
V = np.array([np.cos(a), 0, np.sin(a)])
n = mesh.npanels
s = len(wakes)
pc = S._precompute(mesh)
wc = S.wake_precompute(wakes)
size = np.sqrt(mesh.area)
C = mesh.centroid


def run(rhs_sign, asm_sign, ksign):
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet',
                       kutta_sign=ksign)
    sigma = (mesh.normal @ V)
    rhs = np.zeros(n + s)
    rhs[:n] = rhs_sign * (Brow @ sigma)
    x = np.linalg.solve(A, rhs)
    mu, muw = x[:n], x[n:]
    vel = np.zeros((n, 3))
    for i in range(n):
        p = C[i]
        v = V + asm_sign * (
            (sigma[:, None] * S._source_vel_rows(
                p, pc['qp'], pc['qw'], mesh.area, size,
                self_idx=i, normal=mesh.normal[i])).sum(axis=0)
            + (mu[:, None] * S._loop_vel_rows(
                p, pc['ed'], pc['emask'])).sum(axis=0)
            + (muw[wc['owner'], None] * S._loop_vel_rows(
                p, wc['ed'], wc['emask'])).sum(axis=0))
        # NOTE: self jump is normal-only; removed by tangentialization below.
        vel[i] = v
    vel = vel - (np.einsum('ij,ij->i', vel, mesh.normal))[:, None] * mesh.normal
    cp = 1 - np.einsum('ij,ij->i', vel, vel)
    f = pressure_forces(mesh, cp, V, sref=6.0)
    return f['CL'], f['CDp']


print('(lifting-line AR6 NACA0012 alpha=4: CL~0.33, CDi~0.006)')
for rhs_sign in [+1, -1]:
    for asm_sign in [+1, -1]:
        for ksign in [+1, -1]:
            cl, cdp = run(rhs_sign, asm_sign, ksign)
            print(f'RHS={rhs_sign:+d} ASM={asm_sign:+d} K={ksign:+d} '
                  f'CL={cl:+.3f} CDp={cdp:+.4f}')
