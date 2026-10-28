"""WING3D_BL ablation: viscous-inviscid interaction on the known anchors.

Same campaign configuration as the entropy work (span 1.6, g80x32x32,
ep8, lifting-line imprint, full potential) so the only difference
between a BL0 and a BL1 row is the integral boundary layer.

For each anchor case:
  * run with WING3D_BL=0 (inviscid baseline),
  * run with WING3D_BL=1 (Cp -> integral BL -> delta* transpiration ->
    re-solve, one coupling iteration per nonlinear epoch),
  * grade both against tests3d/validation/known_anchors.json.

Honest-reporting rules: no tuning of alpha, Re or the transition point
to fit the anchors. Transition is FIXED at x/c = 0.05 (WING3D_BL_XTR)
and Re_c = 1.8e7 (WING3D_RE), the latter taken from the OpenFOAM
anchor case freestream (T=300 K, p=100 kPa, U=300 m/s, chord 1 m).

Usage: python tsd_bl.py [smoke|full]
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

ANCHORS = os.path.join(REPO, 'tests3d', 'validation', 'known_anchors.json')
OUT = os.path.join(REPO, 'tests3d', 'validation', 'results')

CASES = [
    ('dfvlr_r4', 0.850),
    ('dfvlr_r4', 0.925),
    ('dfvlr_r4', 1.000),
    ('rae2822', 0.850),
]


def anchor_cl(geo, M):
    """Anchor CL for (geo, M).

    dfvlr_r4 uses fine_125k_dfvlr (the DFVLR fine-branch anchors),
    rae2822 uses the 86k corpus -- same grading as the entropy work.
    Keys are matched as floats, not strings (note "0.925" does not
    round-trip through '%.2f').
    """
    with open(ANCHORS) as f:
        A = json.load(f)
    src = 'fine_125k_dfvlr' if geo == 'dfvlr_r4' else 'corpus_86k'
    tbl = A[src][geo] if src == 'corpus_86k' else A[src]
    for k, v in tbl.items():
        if abs(float(k) - M) < 1e-9:
            return float(v), src
    raise KeyError('%s M%.3f not in %s' % (geo, M, src))


def run_one(geo, M, bl, epochs=8, span=1.6):
    from tsd_campaign import build_profile, build_section
    from tsd_m_sup import run_sup
    os.environ['WING3D_BL'] = '1' if bl else '0'
    prof = build_profile(geo)
    sec = build_section(geo)
    t0 = time.time()
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=False,
                            profile=prof, section=sec)
    dt = time.time() - t0
    cl = float(L['CL'])
    a, src = anchor_cl(geo, M)
    rec = {'geo': geo, 'M': M, 'bl': int(bl), 'CL': cl,
           'anchor': a, 'anchor_src': src,
           'err_pct': 100.0 * (cl - a) / a, 'sec': dt,
           'epochs': epochs, 'span': span,
           'grid': '80x32x32', 'Re': float(os.environ.get('WING3D_RE',
                                                          '1.8e7')),
           'xtr': float(os.environ.get('WING3D_BL_XTR', '0.05'))}
    # centreline section cl + peak suction for the "does it weaken the
    # suction" check
    cl_sec = float(np.trapezoid(L['cpl'][ix, j0] - L['cpu'][ix, j0], xc))
    rec['cl_section'] = cl_sec
    rec['cp_min'] = float(L['cpu'][ix, j0].min())
    rec['x_cp_min'] = float(xc[int(np.argmin(L['cpu'][ix, j0]))])
    if bl and 'bl' in L:
        bl = L['bl']
        rec['bl_dstar_max'] = float(max(np.max(bl['dstar_u']),
                                        np.max(bl['dstar_l'])))
        rec['bl_dstar_u'] = np.asarray(bl['dstar_u']).tolist()
        rec['bl_hist'] = [{'dstar_max': h['dstar_max'],
                           'change': h['change'],
                           'n_sep': h['n_sep']} for h in bl['hist']]
        rec['bl_sep'] = bl['sep']
    return rec


def main(mode='full'):
    epochs = 2 if mode == 'smoke' else 8
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for geo, M in CASES:
        for bl in (0, 1):
            tag = 'bl%d' % bl
            path = os.path.join(OUT, 'bl_%s_M%.3f_%s_ep%d.json' % (
                geo, M, tag, epochs))
            if os.path.exists(path):
                with open(path) as f:
                    rec = json.load(f)
                print('cached %s' % path, flush=True)
            else:
                print('=== %s M=%.3f BL=%d (%d ep) ===' % (geo, M, bl,
                                                           epochs),
                      flush=True)
                rec = run_one(geo, M, bl, epochs=epochs)
                with open(path, 'w') as f:
                    json.dump(rec, f, indent=1)
            rows.append(rec)
            print('RESULT %s M%.3f BL%d CL=%.5f d=%+.2f%% (%.0fs)' % (
                geo, M, bl, rec['CL'], rec['err_pct'], rec['sec']),
                flush=True)
    # summary table, BL0 vs BL1
    print('\n| case | M | anchor | BL0 CL (err) | BL1 CL (err) | delta |')
    print('|---|---|---|---|---|---|')
    for i in range(0, len(rows), 2):
        a, b = rows[i], rows[i + 1]
        print('| %s | %.3f | %.6f | %.5f (%+.2f%%) | %.5f (%+.2f%%) | '
              '%+.2f pp |' % (
                  a['geo'], a['M'], a['anchor'],
                  a['CL'], a['err_pct'], b['CL'], b['err_pct'],
                  b['err_pct'] - a['err_pct']))
    with open(os.path.join(OUT, 'bl_summary_ep%d.json' % epochs), 'w') as f:
        json.dump(rows, f, indent=1)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'full')
