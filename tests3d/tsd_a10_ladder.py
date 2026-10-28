"""A10 ladder: free circulation closed by TANGENCY + Morino wake.

Same harness as tsd_a9_ladder.py, new closure (WING3D_KUTTA=4), and the
NEW anchors (B1 rebuilt OpenFOAM as free-tip AR1.6; the old ones are stale).

  python tsd_a10_ladder.py A     low-Mach lifting-line check (M0.3)
  python tsd_a10_ladder.py B     transonic CL(M) dfvlr 0.85/0.925/1.00
  python tsd_a10_ladder.py C     dfvlr M0.85 with WING3D_BL=1
  python tsd_a10_ladder.py D     legacy free_k side by side

Honesty: only the closure flag, the measured row weight and the physics
gates move between rows; alpha (1.25 deg), Re, transition, grid are fixed.
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

# B1 free-tip AR1.6 anchors (rebuilt 2026-10-28; replaces the 125k-fine
# campaign values 0.3210697 / 0.3629209 / 0.3037588)
ANCHORS = {'0.850': 0.27186, '0.925': 0.26253, '1.000': 0.22694}
ANCHORS_RAE = {'0.850': 0.17143}
OUT = os.path.join(REPO, 'tests3d')


def tag():
    s = 'k%s' % os.environ.get('WING3D_KUTTA', '0')
    for e in ('WING3D_ENTROPY', 'WING3D_BL', 'WING3D_MORINO'):
        if os.environ.get(e, '0') == '1':
            s += '_' + e.split('_')[1].lower()
    return s


def one(M, geo='dfvlr_r4', epochs=8, span=1.6, free_k=False, use_profile=True,
        note='', prefix='a10'):
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
           'morino': os.environ.get('WING3D_MORINO', '1'),
           'alpha_deg': 1.25,
           'free_k': bool(free_k), 'profile': use_profile,
           'CL': float(L['CL']), 'secs': time.time() - t0, 'note': note}
    a = (ANCHORS_RAE if geo == 'rae2822' else ANCHORS).get('%.3f' % M)
    if a:
        rec['anchor'] = a
        rec['err_pct'] = 100.0 * (rec['CL'] - a) / a
    if 'kutta' in L:
        G = np.asarray(L['kutta']['Gamma'])
        y = np.asarray(L['kutta']['y'])
        rec['Gamma'] = G.tolist()
        rec['y'] = y.tolist()
        rec['w'] = L['kutta']['w']
        rec['hist'] = L['kutta']['hist']
        if len(G) > 2:
            gmid = float(G[len(G) // 2])
            rec['gamma_mid'] = gmid
            # ellipticity: tip / mid ratio of |Gamma| (elliptic ~0.30
            # at AR1.6; span-flat would be ~1.0)
            rec['gamma_tip_ratio'] = float(
                np.mean([abs(G[0]), abs(G[-1])]) / max(abs(gmid), 1e-12))
    if 'bl' in L:
        rec['bl_hist'] = [{'dstar_max': h['dstar_max'], 'change': h['change'],
                           'n_sep': h['n_sep']} for h in L['bl']['hist']]
        rec['bl_dstar_max'] = float(max(np.max(L['bl']['dstar_u']),
                                        np.max(L['bl']['dstar_l'])))
    name = '%s_%s_M%.3f_%s_ep%d%s.json' % (
        prefix, geo, M, tag(), epochs, '_noprof' if not use_profile else '')
    path = os.path.join(OUT, name)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    extra = ''
    if 'err_pct' in rec:
        extra = '  err=%+.2f%% vs %.5f' % (rec['err_pct'], rec['anchor'])
    if 'gamma_mid' in rec:
        extra += '  G_mid=%.4f  tip/mid=%.2f' % (
            rec['gamma_mid'], rec['gamma_tip_ratio'])
    print('RESULT %s  M=%.3f CL=%.5f%s  (%s) -> %s' % (
        tag(), M, rec['CL'], extra, '%.0fs' % rec['secs'], name),
        flush=True)
    return rec


def main(ladder):
    geo = os.environ.get('A10_GEO', 'dfvlr_r4')
    ep = int(os.environ.get('A10_EP', '8'))
    if ladder == 'A':
        one(0.3, epochs=int(os.environ.get('A10_EP', '4')),
            use_profile=False, note='low-Mach symmetric / LL check')
        one(0.3, epochs=int(os.environ.get('A10_EP', '4')),
            use_profile=True, note='low-Mach dfvlr section')
    elif ladder == 'B':
        for M in (0.85, 0.925, 1.0):
            one(M, epochs=ep)
    elif ladder == 'C':
        one(0.85, epochs=ep, note='BL on top of free circulation')
    elif ladder == 'D':
        for M in (0.85,):
            one(M, epochs=int(os.environ.get('A10_EP', '4')), free_k=True,
                note='free_k (legacy k-sheet)')
    elif ladder == 'E':
        one(0.85, geo='rae2822', epochs=ep, note='rae2822 new anchor')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'A')
