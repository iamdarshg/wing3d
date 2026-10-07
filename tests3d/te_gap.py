import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
import wing3d.geometry as G
from wing3d.solver import solve, build_wake_from_meta, assemble
from wing3d.forces import pressure_forces

orig = G.naca4_section


def run(gap):
    def blunt(code='0012', n=60, cosine=True):
        loop = orig(code, n, cosine).copy()
        loop[0, 1] += gap / 2
        loop[-1, 1] -= gap / 2
        return loop
    G.build_wing.__globals__['naca4_section'] = blunt
    mesh = G.build_wing('0012', span=6.0, chord=1.0, n_chord=24, n_span=12)
    wakes = build_wake_from_meta(mesh)
    A, Brow = assemble(mesh, wakes, kutta_mode='doublet')
    a = np.radians(4)
    V = np.array([np.cos(a), 0, np.sin(a)])
    res = solve(mesh, wakes, V, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], V, sref=6.0)
    return f['CL'], f['CDp']


for gap in [0.0, 0.001, 0.002, 0.005]:
    cl, cdp = run(gap)
    print(f'te_gap={gap}: CL={cl:+.4f} CDp={cdp:+.5f}')
print('(ref LL 0.33, VLM 0.39)')
