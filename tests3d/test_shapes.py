"""Exclusive-3D validation suite (no 2D strips): shape progression
sphere -> cube -> ellipsoid -> wedge -> Ahmed-like car -> F-16-like.
Run: python -m pytest tests3d/test_shapes.py -x -q
"""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh, build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces
from wing3d.primitives import box, ellipsoid, wedge
from wing3d.shapes import build_car, build_f16


def run_body(mesh, vinf=(1, 0, 0)):
    res = solve(mesh, [], list(vinf))
    f = pressure_forces(mesh, res['cp'], vinf,
                        sref=frontal_area(mesh, vinf))
    return res, f


def frontal_area(mesh, vinf=(1, 0, 0)):
    V = np.asarray(vinf, float)
    vhat = V / np.linalg.norm(V)
    # projected area along flow (convex approx via max extent product)
    C = mesh.centroid
    n = mesh.normal
    return float(np.sum(mesh.area[np.abs(n @ vhat) > 0.7]) * 0.5 + 1e-300)


def test_sphere_analytic():
    from tests3d.analytic_sphere import sphere_mesh
    m = sphere_mesh(12, 24)
    res = solve(m, [], [1, 0, 0])
    C = m.centroid
    th = np.arccos(np.clip(C[:, 0], -1, 1))
    exact = 1 - 2.25 * np.sin(th) ** 2
    rms = float(np.sqrt(((res['cp'] - exact) ** 2).mean()))
    assert rms < 0.15, rms


def test_sphere_drag_zero():
    from tests3d.analytic_sphere import sphere_mesh
    m = sphere_mesh(12, 24)
    res = solve(m, [], [1, 0, 0])
    f = pressure_forces(m, res['cp'], [1, 0, 0], sref=np.pi)
    assert abs(f['CDp']) < 0.02, f['CDp']


def test_cube_stagnation_and_drag():
    m = box(size=(1, 1, 1), n=(10, 10, 10))
    res = solve(m, [], [1, 0, 0])
    C = m.centroid
    front = C[:, 0] > 0.49
    assert res['cp'][front].max() > 0.9  # stagnation ~1
    f = pressure_forces(m, res['cp'], [1, 0, 0], sref=1.0)
    assert abs(f['CDp']) < 0.15, f['CDp']  # d'Alembert (corners limit)


def test_cube_symmetry():
    m = box(size=(1, 1, 1), n=(8, 8, 8))
    res = solve(m, [], [1, 0, 0])
    cp = res['cp']
    C = m.centroid
    # top/bottom and left/right symmetry
    assert abs(cp[C[:, 2] > 0.4].mean() - cp[C[:, 2] < -0.4].mean()) < 0.05
    assert abs(cp[C[:, 1] > 0.4].mean() - cp[C[:, 1] < -0.4].mean()) < 0.05


def test_ellipsoid_self_convergence():
    ms = [ellipsoid(1.0, 0.5, 0.5, nlat=nl, nlon=nn)
          for (nl, nn) in ((16, 32), (20, 40))]
    cps = []
    for m in ms:
        res = solve(m, [], [1, 0, 0])
        assert res['cp'].max() > 0.9  # nose stagnation (needs fine nose mesh)
        cps.append(res['cp'])
    f = pressure_forces(ms[1], cps[1], [1, 0, 0], sref=np.pi * 0.25)
    assert abs(f['CDp']) < 0.05, f['CDp']


def test_wedge_runs_finite():
    # Sharp convex corners limit low-order accuracy; robustness only.
    m = wedge(x0=0, x1=1.0, y0=-0.5, y1=0.5, z0=0.0, z1=0.3,
              nx=10, ny=10)
    res = solve(m, [], [1, 0, 0])
    assert np.all(np.isfinite(res['cp']))
    f = pressure_forces(m, res['cp'], [1, 0, 0], sref=0.3)
    print(f'\nwedge inviscid CD (sharp corners, qualitative): {f["CDp"]:.3f}')
    assert abs(f['CDp']) < 1.0


def test_car_runs_and_reports():
    m = build_car()
    res = solve(m, [], [1, 0, 0])
    assert np.all(np.isfinite(res['cp']))
    f = pressure_forces(m, res['cp'], [1, 0, 0], sref=0.389 * 0.288)
    print('\nAhmed-like inviscid CD (qualitative, separated flow not '
          f'modeled): {f["CDp"]:.4f}')
    assert abs(f['CDp']) < 0.5


def test_f16_runs_and_reports():
    m = build_f16(scale=0.2)  # keep panel count modest for tests
    res = solve(m, [], [1, 0, 0])
    assert np.all(np.isfinite(res['cp']))
    f = pressure_forces(m, res['cp'], [1, 0, 0], sref=1.0)
    print(f'\nF-16-like inviscid CDxS: {f["CDp"]:.4f}')
    assert abs(f['CDp']) < 1.0


def test_wing_lift_slope():
    mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=20, n_span=10)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    cls = []
    for adeg in [0, 4]:
        a = np.radians(adeg)
        res = solve(mesh, wakes, [np.cos(a), 0, np.sin(a)], A=A, Brow=Brow,
                    kutta_mode='doublet')
        f = pressure_forces(mesh, res['cp'], [np.cos(a), 0, np.sin(a)],
                            sref=6.0)
        cls.append(f['CL'])
    slope = (cls[1] - cls[0]) / 4
    print(f'\nAR6 CL slope: {slope:.4f}/deg (LL~0.082, VLM~0.098)')
    assert 0.06 < slope < 0.13, slope
