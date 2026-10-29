"""A16 TASK 3: per-span sections -- rectangular equivalence and a
representative flying-wing planform.

  python a16_planform.py rect    rectangular built from the planform path
                                 vs the stock single-section case
  python a16_planform.py fw      representative flying wing: convergence,
                                 CL, span loading
LOW-RES: g40x32x16 (the box is at <1 GB).
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from a16_common import make_solver, run, span_loading          # noqa
from wing3d.geometry import nonplanar_planform, flying_wing_planform  # noqa

GRID = tuple(int(v) for v in os.environ.get('A16_GRID', '40,32,16').split(','))
EP = int(os.environ.get('A16_EP', '1'))
SPAN = float(os.environ.get('A16_SPAN', '1.6'))
OUT = os.path.join(REPO, 'tests3d', 'a16_planform.json')


def _cl_circ(s):
    y, G = span_loading(s)
    return (2.0 * float(np.sum(np.asarray(G, float)
                               * s.dy[np.asarray(s.wing_jy)])) / s.S_ref(),
            y, np.asarray(G, float))


def _run(tag, plan=None, grid=None, ep=None):
    g = GRID if grid is None else grid
    s = make_solver(g[0], g[1], g[2], span=SPAN, alpha=1.25, M=0.3)
    if plan is not None:
        s.set_planform(plan)
    L = run(s, epochs=EP if ep is None else ep)
    clc, y, G = _cl_circ(s)
    gm = s._wing_geom()
    nact = int(gm['mask'].sum())
    chi = gm['zone'] if 'zone' in gm else None
    print('%-14s grid=%dx%dx%d wi=%2d ng=%2d active=%3d  CL=%.5f  '
          'CL_circ=%.5f' % (tag, g[0], g[1], g[2], len(s.wing_ix),
                            len(s.wing_jy), nact, L['CL'], clc), flush=True)
    return {'tag': tag, 'grid': '%dx%dx%d' % g, 'CL': float(L['CL']),
            'CL_circ': clc, 'n_active': nact, 'S_ref': s.S_ref(),
            'eta': (np.abs(y) / (0.5 * SPAN)).tolist(),
            'Gamma': G.tolist(), 'y': y.tolist(),
            'n_ctrl': int(s.n_ctrl), 'mode': os.environ.get('WING3D_KUTTA')}


def rect():
    stock = _run('stock')
    # same rectangle through the multi-section path: constant chord, no
    # sweep, no twist, no dihedral -> must reproduce the stock case
    plan = nonplanar_planform(
        span=SPAN, n_st=41, chord_fn=lambda e: np.ones_like(e),
        xle_fn=lambda e: np.zeros_like(e), z_fn=lambda e: np.zeros_like(e),
        twist_fn=lambda e: np.zeros_like(e))
    alt = _run('planform-rect', plan)
    dCL = abs(alt['CL'] - stock['CL']) / abs(stock['CL'])
    dCc = abs(alt['CL_circ'] - stock['CL_circ']) / abs(stock['CL_circ'])
    print('  rel|dCL|=%.3e   rel|dCL_circ|=%.3e   (accept <=1e-3)'
          % (dCL, dCc))
    res = {'rect': {'stock': stock, 'planform': alt,
                    'rel_dCL': dCL, 'rel_dCL_circ': dCc,
                    'pass': bool(dCL <= 1e-3 and dCc <= 1e-3)}}
    _save(res)
    return 0 if res['rect']['pass'] else 1


def fw():
    plan = flying_wing_planform(span=SPAN, chord_root=0.62, taper=0.18,
                                eta_crank=0.22, sweep_le_deg=32.0,
                                sweep_le_in_deg=8.0, dihedral_deg=0.0,
                                washout_deg=3.0, n_st=41,
                                t_centre=0.16, t_outer=0.10)
    print('  planform: S=%.4f  AR=%.2f  chord root %.3f tip %.3f'
          % (np.trapezoid(plan['chord'], plan['y']),
             SPAN ** 2 / np.trapezoid(plan['chord'], plan['y']),
             plan['chord'][len(plan['y']) // 2], plan['chord'][0]))
    rows = []
    for ep in (1, 2, 3):
        r = _run('fw ep%d' % ep, plan, ep=ep)
        rows.append(r)
    y = np.asarray(rows[-1]['y'], float)
    G = np.asarray(rows[-1]['Gamma'], float)
    eta = np.abs(y) / (0.5 * SPAN)
    o = np.argsort(eta)
    print('  span loading (eta, Gamma):')
    print('    ' + '  '.join('%.2f:%.4f' % (eta[i], G[i]) for i in o))
    gl = 0.5 * (G[o][:-1] + G[o][1:])
    print('  Gamma(centre)/Gamma(max)=%.3f  tip/centre=%.3f'
          % (G[o][-1] / max(np.abs(G).max(), 1e-12),
             G[o][0] / max(G[o][-1], 1e-12)))
    _save({'flying_wing': {'plan': {k: (v.tolist() if isinstance(
        v, np.ndarray) else v) for k, v in plan.items()
        if k in ('y', 'xle', 'chord', 'z', 'twist_deg', 'note')},
        'runs': rows}})
    return 0


def _save(rec):
    old = {}
    if os.path.exists(OUT):
        try:
            with open(OUT) as f:
                old = json.load(f)
        except Exception:
            old = {}
    old.update(rec)
    with open(OUT, 'w') as f:
        json.dump(old, f, indent=1)
    print('wrote', OUT)


if __name__ == '__main__':
    os.environ.setdefault('WING3D_KUTTA', '6')
    what = sys.argv[1] if len(sys.argv) > 1 else 'rect'
    sys.exit({'rect': rect, 'fw': fw}[what]())
