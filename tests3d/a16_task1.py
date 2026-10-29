"""A16 TASK 1: what does the planar-surface assumption cost?

Three measurements, cheaper-first:

  Z    plane-thickness / snapping sensitivity of the CURRENT model:
       CL at fixed (nx,ny) while the z-node distribution around the body
       plane is varied (nz, the centre-cluster amplitude, and a rigid
       z-offset that puts z=0 between nodes).  Every CL here is a fully
       coupled Kutta solve (mode 6 by default), so the number reported is
       the solver's own answer, not a proxy.

  V    the non-planar reference (tests3d/a16_vlm.py): dihedral 2/5/10 deg,
       spanwise height variation (arch/bump), each as dCL/CL.

  D    the same dihedral geometry handed to the TSD solver (only valid
       once WING3D_NONPLANAR exists; skipped otherwise).

  python a16_task1.py Z|V|all
LOW-RES: g40x24x16 unless A16_GRID says otherwise (the box is at <1 GB).
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from a16_common import make_solver, dzcore, zgrid, run, span_loading  # noqa
import a16_vlm as V                                                  # noqa

GRID = tuple(int(v) for v in os.environ.get('A16_GRID', '40,24,16').split(','))
EP = int(os.environ.get('A16_EP', '1'))
OUT = os.path.join(REPO, 'tests3d', 'a16_task1.json')


def _one(tag, nz=GRID[2], zdist=None, note='', grid=None):
    g = GRID if grid is None else grid
    t0 = time.time()
    s = make_solver(g[0], g[1], nz, zdist=zdist)
    zb, za = dzcore(s)
    L = run(s, epochs=EP)
    y, G = span_loading(s)
    rec = {'tag': tag, 'grid': '%dx%dx%d' % (g[0], g[1], nz),
           'mode': os.environ.get('WING3D_KUTTA', '6'),
           'dz_below': float(zb), 'dz_above': float(za),
           'dz_core': float(min(zb, za)),
           'z0_offset': float(s.z[s.k0]),
           'CL': float(L['CL']), 'CD': float(L['CD']),
           'Gamma_mid': float(G[len(G) // 2]) if len(G) else float('nan'),
           'secs': round(time.time() - t0, 1), 'note': note}
    print('  %-22s %-11s dz=%.4f/%.4f CL=%.5f (%+.2f%%) %.0fs' % (
        tag, rec['grid'], zb, za, rec['CL'], 0.0, rec['secs']), flush=True)
    return rec


def zstudy():
    recs = []
    # 1) nz sweep at fixed nx,ny: changes dz near the wing only
    for nz in (12, 16, 24, 32):
        recs.append(_one('nz=%d' % nz, nz=nz))
    # 2) cluster sweep at fixed nz: dz_core changes 2.5x at the same nz
    for cl in (0.3, 0.75, 0.9):
        recs.append(_one('cluster=%.2f' % cl, nz=16,
                         zdist=zgrid(16, cluster=cl)))
    # 3) rigid z offset: body plane no longer a node (snapping error)
    z0 = zgrid(16)
    for frac in (0.0, 0.25, 0.5):
        d = frac * 0.0406
        recs.append(_one('zoff=%.2fdz' % frac, nz=16, zdist=z0 + d))
    return recs


def vstudy():
    out = {'dihedral': [], 'height': []}
    base = V.vlm(V.stations_rect(span=1.6, ns=96), 1.25, 0.0)
    for g in (0.0, 2.0, 5.0, 10.0, 20.0):
        st = V.stations_rect(span=1.6, ns=96,
                             zfun=lambda yy, g=g: np.abs(yy) * np.tan(
                                 np.radians(g)))
        r = V.vlm(st, 1.25, 0.0)
        out['dihedral'].append({
            'deg': g, 'CL': r['CL'],
            'dCL_pct': 100 * (r['CL'] - base['CL']) / base['CL'],
            'cos': float(np.cos(np.radians(g)))})
    for h in (0.05, 0.10, 0.20):
        st = V.stations_rect(
            span=1.6, ns=96,
            zfun=lambda yy, h=h: h * np.cos(np.pi * yy / 1.6))
        r = V.vlm(st, 1.25, 0.0)
        out['height'].append({'kind': 'arch', 'h': h, 'CL': r['CL'],
                              'dCL_pct': 100 * (r['CL'] - base['CL'])
                              / base['CL']})
    for h in (0.05, 0.10):
        st = V.stations_rect(
            span=1.6, ns=96,
            zfun=lambda yy, h=h: h * (1.0 - (2 * yy / 1.6) ** 2) ** 2)
        r = V.vlm(st, 1.25, 0.0)
        out['height'].append({'kind': 'bump', 'h': h, 'CL': r['CL'],
                              'dCL_pct': 100 * (r['CL'] - base['CL'])
                              / base['CL']})
    out['base_CL'] = base['CL']
    print(json.dumps(out, indent=1))
    return out


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    res = {}
    if what in ('Z', 'all'):
        print('--- Z: plane-thickness / snapping CL sensitivity (%dx%dx* ) '
              '---' % GRID[:2], flush=True)
        res['Z'] = zstudy()
        ref = [r for r in res['Z'] if r['tag'] == 'nz=16'][0]['CL']
        for r in res['Z']:
            print('  %-14s CL=%.5f  dCL=%+.3f%%  dz_core=%.4f' % (
                r['tag'], r['CL'], 100 * (r['CL'] - ref) / ref,
                r['dz_core']))
    if what in ('V', 'all'):
        print('--- V: non-planar VLM reference (relative effects) ---')
        res['V'] = vstudy()
    old = {}
    if os.path.exists(OUT):
        try:
            with open(OUT) as f:
                old = json.load(f)
        except Exception:
            old = {}
    old.update(res)
    with open(OUT, 'w') as f:
        json.dump(old, f, indent=1)
    print('wrote', OUT)
