"""Full-potential driver: exact prescription + line relaxation."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box
from scipy.spatial import cKDTree


def run_fp(Minf, alpha_deg=1.25, nx=80, ny=32, nz=32, itmax=400):
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(alpha_deg)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    meta = mesh.meta
    nl, ns = meta['n_loop'], meta['n_span']
    strips, _ = wing_strips_from_structured(nl, ns)
    info = [(up, lo, meta['te_points'][j], meta['te_points'][j + 1])
            for (up, lo, j) in strips]
    wakes_long = build_wake(mesh, info, length=9.0, n_panels=6)
    g = {'y': np.array([mesh.centroid[j * 48][1] for j in range(24)]),
         'gamma': np.array(res['muw'])}
    s = TSDSolver(Minf=Minf, alpha_deg=alpha_deg, nx=nx, ny=ny, nz=nz,
                  omega=0.8, wake_gamma=g, linear=False, fullpot=True)
    C = mesh.centroid
    mu = np.asarray(res['mu'])
    up_idx = C[:, 2] > 0.001
    lo_idx = C[:, 2] < -0.001
    Tu, Tl = cKDTree(C[up_idx][:, :2]), cKDTree(C[lo_idx][:, :2])
    mu_up, mu_lo = mu[up_idx], mu[lo_idx]
    xc = s.xc[s.wing_ix]
    yc = s.yc[s.wing_jy]
    XX, YY = np.meshgrid(xc, yc, indexing='ij')
    pts = np.stack([XX.ravel(), YY.ravel()], axis=-1)
    s.wing_jump = (mu_up[Tu.query(pts)[1]] -
                   mu_lo[Tl.query(pts)[1]]).reshape(XX.shape)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    s.wake_scale = 1.0
    s.solve(itmax=itmax, tol=1e-7, verbose=True, ramp=100,
            freeze_every=50, init_linear=False, use_lines=True)
    L = s.fp_loads()
    return s, L


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    ad = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    s, L = run_fp(Minf, alpha_deg=ad)
    print('FULLPOT M%.2f a=%.2f: CL=%.4f CD=%.5f' % (Minf, ad, L['CL'],
                                                     L['CD']))
