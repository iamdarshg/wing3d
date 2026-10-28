"""A9 ladder: free bound-vortex circulation closed by Kutta.

Runs the validation ladder of the Kutta-closure work and writes one JSON
per case to tests3d/a9_*.json (never to tests3d/validation/*.json --
that tree is owned elsewhere).

  python tsd_a9_ladder.py A              low-Mach lifting-line check (M0.3)
  python tsd_a9_ladder.py B              dfvlr M0.85/0.925/1.00 vs anchors
  python tsd_a9_ladder.py C              dfvlr M0.85 with WING3D_BL=1
  python tsd_a9_ladder.py D              free_k (old closure) side by side

Env is passed straight through: WING3D_KUTTA=1 (thickness-only + bound
sheet) or 2 (classical free wake), WING3D_ENTROPY, WING3D_BL, ...
Honesty: no alpha/Re/transition/grid tuning; only the formulation flag
and the measured Kutta-row weight move between rows.
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

ANCHORS = {'0.850': 0.3210697, '0.925': 0.3629209, '1.000': 0.3037588}
OUT = os.path.join(REPO, 'tests3d')


def tag():
    m = os.environ.get('WING3D_KUTTA', '0')
    s = 'k%s' % m
    for e in ('WING3D_ENTROPY', 'WING3D_BL'):
        if os.environ.get(e, '0') == '1':
            s += '_' + e.split('_')[1].lower()
    return s


def one(M, geo='dfvlr_r4', epochs=8, span=1.6, free_k=False, use_profile=True,
        note=''):
    from tsd_campaign import build_profile, build_section
    from tsd_coupled import build_case
    from tsd_m_sup import run_sup
    prof = build_profile(geo) if use_profile else None
    sec = build_section(geo)
    t0 = time.time()
    base, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              span, section=sec, profile=prof)
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=free_k,
                            profile=prof, section=sec, base=base)
    rec = {'geo': geo, 'M': M, 'grid': '80x32x32', 'epochs': epochs,
           'span': span, 'mode': os.environ.get('WING3D_KUTTA', '0'),
           'entropy': os.environ.get('WING3D_ENTROPY', '0'),
           'bl': os.environ.get('WING3D_BL', '0'),
           'free_k': bool(free_k), 'profile': use_profile,
           'CL': float(L['CL']), 'secs': time.time() - t0, 'note': note}
    a = ANCHORS.get('%.3f' % M)
    if a:
        rec['anchor'] = a
        rec['err_pct'] = 100.0 * (rec['CL'] - a) / a
    if 'kutta' in L:
        rec['Gamma'] = np.asarray(L['kutta']['Gamma']).tolist()
        rec['y'] = np.asarray(L['kutta']['y']).tolist()
        rec['w'] = L['kutta']['w']
        rec['hist'] = L['kutta']['hist']
    if 'bl' in L:
        rec['bl_hist'] = [{'dstar_max': h['dstar_max'], 'change': h['change'],
                           'n_sep': h['n_sep']} for h in L['bl']['hist']]
        rec['bl_dstar_max'] = float(max(np.max(L['bl']['dstar_u']),
                                        np.max(L['bl']['dstar_l'])))
    name = 'a9_%s_M%.3f_%s_ep%d%s.json' % (
        geo, M, tag(), epochs, '_noprof' if not use_profile else '')
    path = os.path.join(OUT, name)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    print('RESULT %s  M=%.3f CL=%.5f%s  (%s) -> %s' % (
        tag(), M, rec['CL'],
        '  err=%+.2f%%' % rec['err_pct'] if 'err_pct' in rec else '',
        '%.0fs' % rec['secs'], name), flush=True)
    return rec


def main(ladder):
    geo = os.environ.get('A9_GEO', 'dfvlr_r4')
    ep = int(os.environ.get('A9_EP', '8'))
    if ladder == 'A':
        # lifting-line check: 4 epochs, span 1.6. Reference (LL imprint,
        # L/D-free): CL 0.0545 at AR1.6 for the symmetric section.
        rows = []
        rows.append(one(0.3, epochs=int(os.environ.get('A9_EP', '4')),
                        use_profile=False, note='low-Mach symmetric / LL check'))
        rows.append(one(0.3, epochs=int(os.environ.get('A9_EP', '4')),
                        use_profile=True, note='low-Mach dfvlr section'))
        for r in rows:
            print(json.dumps({k: r.get(k) for k in
                              ('M', 'CL', 'profile', 'mode', 'secs')}))
    elif ladder == 'B':
        for M in (0.85, 0.925, 1.0):
            one(M, epochs=ep)
    elif ladder == 'C':
        one(0.85, epochs=ep, note='BL on top of free circulation')
    elif ladder == 'D':
        for M in (0.85, 1.0):
            one(M, epochs=ep, free_k=True, note='free_k (legacy k-sheet)')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'A')
