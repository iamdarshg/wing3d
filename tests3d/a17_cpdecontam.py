"""A17 acceptance: does WING3D_CPDECONTAM collapse the dz floor?

A  dz sweep  nz = 12/16/24/32, M0.3 AR1.6 g40x24x16 stock, KUTTA=6,
   OFF and ON.  OFF must reproduce A16's recorded
   0.07029/0.05541/0.05358/0.04380 bit-identically.
B  rigid z-offset 0.25 dz / 0.5 dz at nz=16 (the body plane stops being a
   node), OFF and ON.
C  campaign regression: g80x32x32, dfvlr_r4, M0.85, ep8,
   WING3D_IMPRINT=ll + WING3D_ENTROPY=1 + KUTTA=0 (the validated config of
   tests3d/sweep_multi.py), OFF vs ON, against the OF anchor 0.32682.

    python a17_cpdecontam.py A|B|C|all
LOW-RES diagnostic grid g40x24x16 for A/B (the box is at ~1.8 GB free).
"""
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

OUT = os.path.join(REPO, 'tests3d', 'a17_cpdecontam.json')
A16 = {12: 0.07029, 16: 0.05541, 24: 0.05358, 32: 0.04380}
ANCHOR = 0.32682


def _solve(nz, zdist=None, M=0.3, grid=(40, 24, 16), ep=1, kutta='6'):
    os.environ['WING3D_KUTTA'] = kutta
    from a16_common import make_solver, run, dzcore
    s = make_solver(grid[0], grid[1], nz, zdist=zdist, M=M)
    L = run(s, epochs=ep, M=M, verbose=False)
    zb, za = dzcore(s)
    return s, L, zb, za


def _off_on(run_kw, tag, phi_hash=True):
    """Same solve twice: CPDECONTAM off then on."""
    rec = {'tag': tag}
    for state in ('off', 'on'):
        os.environ['WING3D_CPDECONTAM'] = '1' if state == 'on' else '0'
        s, L, zb, za = _solve(**run_kw)
        rec['CL_' + state] = float(L['CL'])
        rec['CD_' + state] = float(L['CD'])
        rec['dz'] = [float(zb), float(za)]
        rec['z0'] = float(s.z[s.k0])
        rec['k0'] = int(s.k0)
        if state == 'off':
            rec['phi_off'] = s.phi.copy()
        elif phi_hash:
            rec['dphi_max'] = float(np.abs(s.phi - rec['phi_off']).max())
        del s
    rec.pop('phi_off', None)
    rec['dCL'] = rec['CL_on'] - rec['CL_off']
    print('  %-12s OFF=%.5f ON=%.5f  dCL=%+.5f  dz=%s z0=%+.4f%s'
          % (tag, rec['CL_off'], rec['CL_on'], rec['dCL'],
             '/'.join('%.4f' % v for v in rec['dz']), rec['z0'],
             '' if 'dphi_max' not in rec else
             '  max|dphi(ON-OFF)|=%.1e' % rec['dphi_max']), flush=True)
    return rec


def study_A():
    print('--- A: dz sweep, g%dx%dx* M0.3 AR1.6 KUTTA=6 ---' % (40, 24),
          flush=True)
    recs = [_off_on(dict(nz=nz), 'nz=%d' % nz) for nz in (12, 16, 24, 32)]
    for r in recs:
        nz = int(r['tag'].split('=')[1])
        r['A16_off'] = A16[nz]
        r['bitident'] = abs(r['CL_off'] - A16[nz]) < 5e-5
    return recs


def study_B():
    print('--- B: rigid z-offset, nz=16 ---', flush=True)
    from a16_common import zgrid
    z0 = zgrid(16)
    recs = [_off_on(dict(nz=16), 'zoff=0.00dz')]
    for frac in (0.25, 0.5):
        recs.append(_off_on(dict(nz=16, zdist=z0 + frac * 0.0406),
                            'zoff=%.2fdz' % frac))
    base = recs[0]
    for r in recs:
        r['d_pct_off'] = 100 * (r['CL_off'] - base['CL_off']) / base['CL_off']
        r['d_pct_on'] = 100 * (r['CL_on'] - base['CL_on']) / base['CL_on']
    for r in recs:
        print('  %-12s zoff dCL: OFF=%+.1f%% ON=%+.1f%%'
              % (r['tag'], r['d_pct_off'], r['d_pct_on']), flush=True)
    return recs


def study_C():
    print('--- C: campaign regression dfvlr_r4 M0.85 g80x32x32 ep8 ---',
          flush=True)
    from sweep_multi import apply_validated_env
    apply_validated_env()
    recs = []
    for state in ('off', 'on'):
        os.environ['WING3D_CPDECONTAM'] = '1' if state == 'on' else '0'
        from tsd_coupled import build_case
        from tsd_m_sup import run_sup
        from sweep_multi import build_profile, build_section, GRID, SPAN, ALPHA
        geo = 'dfvlr_r4'
        prof, sec = build_profile(geo), build_section(geo)
        s, _r, _m = build_case(0.3, ALPHA, False, True, None, GRID[0], GRID[1], GRID[2],
                               SPAN, section=sec, profile=prof)
        L, _xc, _ix, _j0 = run_sup(0.85, span=SPAN, epochs=8, free_k=False,
                                   profile=prof, section=sec, base=s,
                                   alpha_deg=ALPHA)
        recs.append({'state': state, 'CL': float(L['CL']),
                     'CD': float(L['CD']),
                     'CL_vs_anchor_pct': 100 * (L['CL'] - ANCHOR) / ANCHOR})
        print('  CPDECONTAM=%s  CL=%.5f  CD=%.5f  vs anchor %.5f: %+.1f%%'
              % (state, L['CL'], L['CD'], ANCHOR,
                 recs[-1]['CL_vs_anchor_pct']), flush=True)
        del s
    return recs


if __name__ == '__main__':
    what = (sys.argv[1] if len(sys.argv) > 1 else 'all').upper()
    res = {}
    if os.path.exists(OUT):
        try:
            res = json.load(open(OUT))
        except Exception:
            res = {}
    for key, fn in (('A', study_A), ('B', study_B), ('C', study_C)):
        if what in (key, 'ALL'):
            res[key] = fn()
    json.dump(res, open(OUT, 'w'), indent=1)
    if 'A' in res:
        off = [r['CL_off'] for r in res['A']]
        on = [r['CL_on'] for r in res['A']]
        print('\nA dz spread:  OFF %.1f%%  (%.5f..%.5f)   ON %.1f%%  '
              '(%.5f..%.5f)' % (
                  100 * (max(off) - min(off)) / np.mean(off), min(off),
                  max(off), 100 * (max(on) - min(on)) / np.mean(on),
                  min(on), max(on)))
        print('   reduction %.1fx   bit-identical OFF: %s' % (
            (max(off) - min(off)) / max(max(on) - min(on), 1e-12),
            all(r['bitident'] for r in res['A'])))
    print('wrote', OUT)
