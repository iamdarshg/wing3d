"""A16 non-planar path: bit-identity gate and dihedral verification.

  python a16_nonplanar.py bit        zsurf = 0 must be BIT-IDENTICAL
  python a16_nonplanar.py dih        uniform dihedral 2/5/10 deg vs the
                                     VLM reference (tests3d/a16_vlm.py)
  python a16_nonplanar.py zs         non-uniform (arch/bump) height
LOW-RES by default (the box is at <1 GB): A16_GRID=40,24,16.
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from a16_common import make_solver, dzcore, run, span_loading      # noqa
import a16_vlm as V                                                # noqa

GRID = tuple(int(v) for v in os.environ.get('A16_GRID', '40,24,16').split(','))
EP = int(os.environ.get('A16_EP', '1'))
SPAN = float(os.environ.get('A16_SPAN', '1.6'))
ALPHA = float(os.environ.get('A16_ALPHA', '1.25'))
MACH = float(os.environ.get('A16_M', '0.3'))
OUT = os.path.join(REPO, 'tests3d', 'a16_nonplanar.json')


def _solve(zsurf=None, grid=None, ep=None, span=SPAN, alpha=ALPHA,
           M=MACH, mode=None, extra=None):
    g = GRID if grid is None else grid
    if mode is not None:
        os.environ['WING3D_KUTTA'] = str(mode)
    s = make_solver(g[0], g[1], g[2], span=span, alpha=alpha, M=M)
    if zsurf is not None:
        s.set_surface(zsurf)
    if extra is not None:
        extra(s)
    t0 = time.time()
    L = run(s, M=M, epochs=EP if ep is None else ep, alpha=alpha, span=span)
    y, G = span_loading(s)
    return {'CL': float(L['CL']), 'CD': float(L['CD']),
            'phi': s.phi.copy(), 'mu': s.mu_b.copy(),
            'y': y, 'G': np.asarray(G, float),
            'secs': round(time.time() - t0, 1),
            'mode': os.environ.get('WING3D_KUTTA', '0'),
            'grid': '%dx%dx%d' % g}


def bit():
    a = _solve(None)
    b = _solve(np.zeros((GRID[0], GRID[1])))
    dCL = abs(b['CL'] - a['CL'])
    dphi = float(np.abs(b['phi'] - a['phi']).max())
    dmu = float(np.abs(b['mu'] - a['mu']).max()) if a['mu'].size else 0.0
    print('planar   CL=%.17g' % a['CL'])
    print('zsurf=0  CL=%.17g' % b['CL'])
    print('max|dCL|=%.3e  max|dphi|=%.3e  max|dmu|=%.3e  (%s, %s)'
          % (dCL, dphi, dmu, a['grid'], a['mode']))
    ok = (dCL == 0.0) and (dphi == 0.0) and (dmu == 0.0)
    print('BIT-IDENTICAL:', ok)
    rec = {'bit': {'CL_planar': a['CL'], 'CL_zsurf0': b['CL'],
                   'dCL': dCL, 'dphi': dphi, 'dmu': dmu,
                   'bit_identical': bool(ok), 'grid': a['grid'],
                   'mode': a['mode'], 'secs': a['secs']}}
    _save(rec)
    return 0 if ok else 1


def _wing_zsurf(g, kind, deg=None, h=None):
    """z_s(x,y) on the (nx,ny) node grid for the test surface shapes."""
    x, y = np.asarray(g[0], float), np.asarray(g[1], float)
    XX, YY = np.meshgrid(x, y, indexing='ij')
    if kind == 'dih':
        return np.abs(YY) * np.tan(np.radians(deg))
    if kind == 'arch':
        return h * np.cos(np.pi * YY / SPAN)
    if kind == 'bump':
        return h * (1.0 - (2 * YY / SPAN) ** 2) ** 2
    if kind == 'sweep_h':      # spanwise LINEAR height ramp (not dihedral:
        return h * (YY / (0.5 * SPAN))    #  port low, starboard high)
    raise ValueError(kind)


def _dyy(y):
    """Span-cell widths for the span stations (edge cells use the end
    difference), the same convention `s.dy[s.wing_jy]` uses."""
    y = np.asarray(y, float)
    if len(y) < 2:
        return np.zeros(len(y))
    d = np.empty(len(y))
    d[1:-1] = 0.5 * (y[2:] - y[:-2])
    d[0] = y[1] - y[0]
    d[-1] = y[-1] - y[-2]
    return d


def dih():
    z0 = _solve(None)
    S = SPAN * 1.0
    z0c = 2.0 * float(np.sum(np.asarray(z0['G'], float)
                             * _dyy(z0['y']))) / S
    rec = {'dihedral': {'flat_CL': z0['CL'], 'flat_CL_circ': z0c,
                        'base': 'planar'}}
    rows = []
    for deg in (2.0, 5.0, 10.0, 20.0):
        s = make_solver(GRID[0], GRID[1], GRID[2], span=SPAN, alpha=ALPHA,
                        M=MACH)
        zs = _wing_zsurf((s.x, s.y), 'dih', deg=deg)
        s.set_surface(zs)
        L = run(s, M=MACH, epochs=EP, alpha=ALPHA, span=SPAN)
        st = V.vlm(V.stations_rect(span=SPAN, ns=96,
                                   zfun=lambda yy, g=deg: np.abs(yy)
                                   * np.tan(np.radians(g))), ALPHA, MACH)
        sf = V.vlm(V.stations_rect(span=SPAN, ns=96), ALPHA, MACH)
        y, G = span_loading(s)
        S = s.S_ref()
        cl_circ = 2.0 * float(np.sum(np.asarray(G, float)
                                     * s.dy[np.asarray(s.wing_jy)])) / S
        row = {'deg': deg, 'CL': float(L['CL']),
               'CL_circ': cl_circ,
               'dCL_pct': 100 * (L['CL'] - z0['CL']) / z0['CL'],
               'dCLc_pct': 100 * (cl_circ - z0c) / z0c,
               'vlm_dCL_pct': 100 * (st['CL'] - sf['CL']) / sf['CL'],
               'cos_pct': 100 * (np.cos(np.radians(deg)) - 1.0),
               'zs_max': float(np.abs(zs).max()),
               'dz': float(min(dzcore(s))),
               'Gmid': float(G[len(G) // 2]) if len(G) else float('nan'),
               'conv': 'ok'}
        rows.append(row)
        print('  dih %4.1f  CL=%.5f (%+.2f%%)  CL_circ=%.5f (%+.3f%%)  '
              'VLM=%+.3f%%  cos=%+.3f%%  zs_max=%.4f (%.2f dz)' % (
                  deg, row['CL'], row['dCL_pct'], cl_circ, row['dCLc_pct'],
                  row['vlm_dCL_pct'], row['cos_pct'], row['zs_max'],
                  row['zs_max'] / row['dz']), flush=True)
    rec['dihedral']['rows'] = rows
    _save(rec)
    return 0


def zs():
    z0 = _solve(None)
    rows = []
    for kind, val in (('arch', 0.1), ('bump', 0.1), ('sweep_h', 0.05),
                      ('sweep_h', 0.10)):
        s = make_solver(GRID[0], GRID[1], GRID[2], span=SPAN, alpha=ALPHA,
                        M=MACH)
        z = _wing_zsurf((s.x, s.y), kind, h=val)
        s.set_surface(z)
        L = run(s, M=MACH, epochs=EP, alpha=ALPHA, span=SPAN)
        if kind == 'sweep_h':
            vz = V.vlm(V.stations_rect(
                span=SPAN, ns=96,
                zfun=lambda yy, h=val: h * (yy / (0.5 * SPAN))), ALPHA,
                MACH)
        elif kind == 'arch':
            vz = V.vlm(V.stations_rect(
                span=SPAN, ns=96,
                zfun=lambda yy, h=val: h * np.cos(np.pi * yy / SPAN)),
                ALPHA, MACH)
        else:
            vz = V.vlm(V.stations_rect(
                span=SPAN, ns=96,
                zfun=lambda yy, h=val: h * (1 - (2 * yy / SPAN) ** 2) ** 2),
                ALPHA, MACH)
        vs = V.vlm(V.stations_rect(span=SPAN, ns=96), ALPHA, MACH)
        rows.append({'kind': kind, 'h': val, 'CL': float(L['CL']),
                     'dCL_pct': 100 * (L['CL'] - z0['CL']) / z0['CL'],
                     'vlm_dCL_pct': 100 * (vz['CL'] - vs['CL']) / vs['CL']})
        print('  %-8s h=%.2f  CL=%.5f  dCL=%+.3f%%  VLM=%+.3f%%' % (
            kind, val, L['CL'], rows[-1]['dCL_pct'],
            rows[-1]['vlm_dCL_pct']), flush=True)
    _save({'zsurf_study': {'flat_CL': z0['CL'], 'rows': rows}})
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
    what = sys.argv[1] if len(sys.argv) > 1 else 'bit'
    sys.exit({'bit': bit, 'dih': dih, 'zs': zs}[what]())
