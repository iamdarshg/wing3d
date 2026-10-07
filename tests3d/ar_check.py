import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import build_wing
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

for span in [6.0, 20.0]:
    ns = 16 if span == 6.0 else 28
    mesh = build_wing('0012', span=span, chord=1.0, n_chord=24, n_span=ns)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    area = span * 1.0
    f = pressure_forces(mesh, res['cp'], V, sref=area)
    nl = mesh.meta['n_loop']
    npl = nl - 1
    dy = span / ns
    jmid = ns // 2
    ids = [jmid * npl + i for i in range(npl)]
    Fz = (-res['cp'][ids, None] * mesh.normal[ids]
          * mesh.area[ids, None]).sum(axis=0)[2]
    # sectional cl = L'/(q*c), L' = (F_raw/dy)*q with F_raw w/o q:
    # cl = F_raw/(dy) since q cancels: cl = Fz_raw/(0.5*1*dy)*0.5... careful:
    # Fz_raw = sum(-Cp n A) [units q*L^2]; L' = Fz_raw*q/dy; cl = L'/(q*1) = Fz_raw/dy
    print(f'span={span} AR={span}: total CL={f["CL"]:.4f} '
          f'midspan sec_cl={Fz / dy:.4f}')
print('(2D cl=0.46, LL AR6=0.33, AR20=0.42)')
