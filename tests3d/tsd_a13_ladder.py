"""A13 ladder: mode 7 = HYBRID tangency row (field-induced w minus the
estimator's analytic bias, plus the A12 kernel), so the circulation is BOTH
physically shaped AND Mach-coupled.

  python tsd_a13_ladder.py A     M0.3 symmetric   (test A: must not regress
                                 A12's 0.0460, and d(row)/d(phi) must be
                                 clearly nonzero where mode 6's was 0)
  python tsd_a13_ladder.py Ac    M0.3 cambered dfvlr (test D)
  python tsd_a13_ladder.py B     transonic dfvlr 0.85/0.925/1.00 (test B)
  python tsd_a13_ladder.py C     dfvlr M0.85 + WING3D_BL=1 (test C)

Same harness as tsd_a12_ladder.py (`one` imported from tsd_a11_ladder) so the
grids, alpha, Re and the anchors are identical; only the closure flag moves.
`one7` is a copy of that recorder that also stores the mode-7 row
decomposition (w_field, its analytic equivalent, the bias B.mu and the
phi/mu coupling report), which the flat A11 record would drop.
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_a11_ladder import ANCHORS, OUT, tag, one    # noqa: E402,F401
from tsd_coupled import build_case                   # noqa: E402
from tsd_m_sup import run_sup                        # noqa: E402
from tsd_campaign import build_profile, build_section  # noqa: E402


def one7(M, geo='dfvlr_r4', epochs=4, span=1.6, use_profile=True,
         note='', prefix='a13'):
    """The A11 recorder + the mode-7 row decomposition."""
    prof = build_profile(geo) if use_profile else None
    sec = build_section(geo)
    t0 = time.time()
    base, _r, _m = build_case(0.3, 1.25, False, True, None, 80, 32, 32,
                              span, section=sec, profile=prof)
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=False,
                            profile=prof, section=sec, base=base)
    rec = {'geo': geo, 'M': M, 'grid': '80x32x32', 'epochs': epochs,
           'span': span, 'mode': os.environ.get('WING3D_KUTTA', '0'),
           'entropy': os.environ.get('WING3D_ENTROPY', '0'),
           'bl': os.environ.get('WING3D_BL', '0'),
           'nc': os.environ.get('WING3D_KUTTA_NC', 'all'),
           'te_kutta': os.environ.get('WING3D_KUTTA_TE', '1'),
           'alpha_deg': 1.25, 'free_k': False, 'profile': use_profile,
           'CL': float(L['CL']), 'CD': float(L['CD']),
           'secs': time.time() - t0, 'note': note}
    a = ANCHORS.get('%.3f' % M)
    if a:
        rec['anchor'] = a
        rec['err_pct'] = 100.0 * (rec['CL'] - a) / a
    K = L['kutta']
    rec['w'] = K.get('w')
    rec['hist'] = K['hist']
    y = np.asarray(K['y'])
    j0 = int(np.argmin(np.abs(y)))
    mu = np.asarray(K['mu'])
    Gsec = np.asarray(K['Gamma_mu'])
    load = np.asarray(K['loading'])
    rec['x_ctrl'] = np.asarray(K['x_ctrl']).tolist()
    rec['mu_mid'] = mu[:, j0].tolist()
    rec['loading_mid'] = load[:, j0].tolist()
    rec['Gamma_mu'] = Gsec.tolist()
    rec['Gamma_mu_mid'] = float(Gsec[j0])
    rec['tang_rows_mid'] = np.asarray(K['tang_rows'])[:, j0].tolist()
    rec['loading_TE_ratio'] = float(
        abs(load[-1, j0]) / max(abs(load[:, j0]).max(), 1e-12))
    G = np.asarray(K['Gamma'])
    rec['y'] = y.tolist()
    rec['Gamma'] = G.tolist()
    if len(G) > 2:
        gmid = float(G[len(G) // 2])
        rec['gamma_mid'] = gmid
        rec['gamma_tip_ratio'] = float(
            np.mean([abs(G[0]), abs(G[-1])]) / max(abs(gmid), 1e-12))
    if 'w_field' in K:
        # A13: the row decomposition. w_field = the estimator on the FIELD,
        # w_eq = the same estimator on the analytic sheet of the SAME mu, so
        # w_field - w_eq is what the field carries beyond the ideal sheet
        # and bias = B.mu is what mode 7 removes.
        wf = np.asarray(K['w_field'])
        we = np.asarray(K['w_eq'])
        bi = np.asarray(K['bias'])
        rec['w_field'] = wf.tolist()
        rec['w_eq'] = we.tolist()
        rec['bias'] = bi.tolist()
        rec['resid_phys'] = (wf - we).tolist()
        rec['coupling'] = K['coupling']
        al = np.radians(rec['alpha_deg'])
        rec['w_field_mid_alpha'] = (wf[:, j0] / al).tolist()
        rec['w_eq_mid_alpha'] = (we[:, j0] / al).tolist()
        rec['bias_mid_alpha'] = (bi[:, j0] / al).tolist()
        rec['resid_phys_mid_alpha'] = ((wf[:, j0] - we[:, j0]) / al).tolist()
    if 'bl' in L:
        rec['bl_hist'] = [{'dstar_max': h['dstar_max'],
                           'change': h['change'],
                           'n_sep': h['n_sep']} for h in L['bl']['hist']]
        rec['bl_dstar_max'] = float(max(np.max(L['bl']['dstar_u']),
                                        np.max(L['bl']['dstar_l'])))
    name = '%s_%s_M%.3f_%s_ep%d%s.json' % (
        prefix, geo, M, tag(), epochs, '' if use_profile else '_noprof')
    path = os.path.join(OUT, name)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    extra = ''
    if 'err_pct' in rec:
        extra = '  err=%+.2f%% vs %.5f' % (rec['err_pct'], rec['anchor'])
    if 'gamma_mid' in rec:
        extra += '  G_mid=%.4f  tip/mid=%.2f' % (
            rec['gamma_mid'], rec['gamma_tip_ratio'])
    if 'coupling' in rec:
        c = rec['coupling']
        extra += '  dphi=%.1f dmu=%.1f mode6_dphi=%.0e' % (
            c['dphi_row'], c['dmu_row'], c['dphi_row_mode6'])
    print('RESULT %s  M=%.3f CL=%.5f%s  (%s) -> %s' % (
        tag(), M, rec['CL'], extra, '%.0fs' % rec['secs'], name), flush=True)
    return rec


def main(ladder):
    ep = int(os.environ.get('A13_EP', os.environ.get('A12_EP', '2')))
    epa = int(os.environ.get('A13_EPA', '4'))
    if ladder == 'A':
        one7(0.3, epochs=epa, use_profile=False, note='mode7 M0.3 symmetric')
    elif ladder == 'Ac':
        one7(0.3, epochs=epa, note='mode7 M0.3 cambered')
    elif ladder == 'B':
        for M in (0.85, 0.925, 1.0):
            one7(M, epochs=ep, note='mode7 transonic')
    elif ladder == 'C':
        one7(0.85, epochs=ep, note='mode7 + BL')
    elif ladder == 'A6':
        # the A12 control through the SAME recorder (no mode-7 keys)
        one(0.3, epochs=epa, use_profile=False, prefix='a13',
            note='mode6 control M0.3 symmetric')
    else:
        raise SystemExit('unknown ladder %s' % ladder)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'A')
