"""TSD (nonlinear) vs OpenFOAM anchor at matched alpha/Mach."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import (solve, build_wake_from_meta, assemble, build_wake,
                           wing_strips_from_structured)
from wing3d.tsd import TSDSolver, tsd_farfield_box


def run_anchor(Minf, alpha_deg=1.25, nx=80, ny=32, nz=32, itmax=1200):
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
                  omega=0.8, wake_gamma=g)
    # panel-mu-exact bound sheet (sign calibrated: mu_up - mu_lo)
    from scipy.spatial import cKDTree
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
    iu = Tu.query(pts)[1]
    il = Tl.query(pts)[1]
    s.wing_jump = (mu_up[iu] - mu_lo[il]).reshape(XX.shape)
    s.farfield = tsd_farfield_box(s, wakes_long, np.array(res['muw']))
    s.wake_scale = 0.0
    s.solve(itmax=itmax, tol=1e-7, verbose=True, ramp=200,
            freeze_every=150, init_linear=True)
    L = s.loads()
    print('TSD M%.2f a=%.2f: CL=%.4f CD=%.5f' % (Minf, alpha_deg,
                                                 L['CL'], L['CD']),
          flush=True)
    # midspan Cp for shock comparison
    jmid = int(np.argmin(np.abs(s.yc)))
    xs, cpu, cpl = [], [], []
    for i in s.wing_ix:
        xs.append(s.xc[i])
        cpu.append(L['cpu'][i, jmid])
        cpl.append(L['cpl'][i, jmid])
    return s, L, np.array(xs), np.array(cpu), np.array(cpl)


if __name__ == '__main__':
    Minf = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8
    s, L, xs, cpu, cpl = run_anchor(Minf)
    print('midspan upper Cp(x):')
    for x, c in zip(xs[::4], cpu[::4]):
        print('  x=%.3f Cp=%.3f' % (x, c))
